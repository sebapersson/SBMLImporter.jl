"""
    load_SBML(path_SBML::AbstractString;
              ifelse_to_callback::Bool=true,
              inline_assignment_rules::Bool=false,
              write_to_file::Bool=false,
              model_as_string::Bool=false,
              check_massaction=true,
              mass_action::Bool=false)

Parse an SBML model into a `ParsedReactionNetwork` and convert SBML events/piecewise to callbacks.

For information on simulating the `ParsedReactionNetwork`, as a `JumpProblem`, a `SDEProblem`, or
an `ODEProblem` see the documentation.

`path_SBML` can be the model as a string if `model_as_string=true`.

## Arguments
- `path_SBML`: File path to a valid SBML file (level 2 or higher).
- `ifelse_to_callback=true`: Whether to rewrite `ifelse` (piecewise) expressions to callbacks; recommended
    for performance.
- `inline_assignment_rules=true`: Whether to inline assignment rules into model equations. Recommended for
    model import speed, however, it will not be possible to access the rule-variable via `sol[:var]`.
- `write_to_file=false`: Whether to write the parsed SBML model to a Julia file in the same directory as the
    SBML file.
- `model_as_string=false`: Whether or not the model (`path_SBML`) is provided as a string.
- `check_massaction=true`: Whether to check if and rewrite a reaction to mass action if possible.
- `mass_action=false`: Determines if reactions are treated as mass-action reactions. This option should
    be set to true **only** if it is certain that the model has mass action reactions, for example
    if the model is an SBML file produced by BioNetGen.

## Returns
- `parsed_rn`: A `ParsedReactionNetwork` struct that can be converted into a `JumpProblem`, a `SDEProblem`, or
    an `ODEProblem`
- `cbset`: Callbackset (events, piecewise etc...) for the model.

## Examples
```julia
# Import and simulate model as a JumpProblem
using SBMLImporter, JumpProcesses
prnbng, cb = load_SBML(path_SBML)
tspan = (0.0, 10.0)
dprob = DiscreteProblem(prnbng.rn, prnbng.u₀, tspan, prnbng.p)
jprob = JumpProblem(prnbng.rn, dprob, Direct())
sol = solve(jprob, SSAStepper(), callback=cb)
```
```julia
# Import and simulate model as a SDE
using SBMLImporter, StochasticDiffEq
prnbng, cb = load_SBML(path_SBML)
tspan = (0.0, 10.0)
sprob = SDEProblem(prnbng.rn, prnbng.u₀, tspan, prnbng.p)
sol = solve(sprob, LambaEM(), callback=cb)
```
```julia
# Import and simulate model as an ODE
using SBMLImporter, ModelingToolkit, OrdinaryDiffEq
prnbng, cb = load_SBML(path_SBML)
sys = convert(ODESystem, prnbng.rn)
oprob = ODEProblem(sys, prnbng.u₀, tspan, prnbng.p, jac=true)
sol = solve(oprob, Rodas5P(), callback=cb)
```
"""
function load_SBML(path_SBML::AbstractString; ifelse_to_callback::Bool = true,
                   inline_assignment_rules::Bool = true, write_to_file::Bool = false,
                   model_as_string::Bool = false, check_massaction::Bool = true,
                   mass_action::Bool = false)
    model_SBML = build_SBML_model(path_SBML; ifelse_to_callback = ifelse_to_callback,
                                  model_as_string = model_as_string,
                                  inline_assignment_rules = inline_assignment_rules,
                                  mass_action = mass_action)
    # If model is written to file save it in the same directory as the SBML-file. Only
    # save if model is not provided as a string (as then there is not file)
    dir_save = _get_dir_save(write_to_file, model_as_string)

    model_SBML_sys = _to_system_syntax(model_SBML, inline_assignment_rules; check_massaction = check_massaction)
    rn, specie_map, parameter_map = _get_reaction_system(model_SBML_sys, model_SBML.name)

    # Build callback functions
    cbset, callback_str = create_callbacks_SBML(rn, model_SBML, model_SBML.name)

    # if model is written to file write the callback
    #= TODO: Fix later
    if write_to_file == true
        path_save = joinpath(dir_save, model_SBML.name * "_callbacks.jl")
        io = open(path_save, "w")
        write(io, callback_str)
        close(io)
        _ = reactionsystem_to_string(parsed_model_SBML, write_to_file, path_save_model,
                                     model_SBML)
    end
    =#

    parsed_rn = ParsedReactionNetwork(rn, specie_map, parameter_map, nothing, nothing)
    return parsed_rn, cbset
end

function _get_dir_save(write_to_file::Bool, model_as_string::Bool)::Union{Nothing, String}
    if write_to_file == true
        dir_save = if model_as_string
            nothing
        else
            joinpath(splitdir(path_SBML)[1], "SBML")
        end
    else
        dir_save = nothing
    end
    if !(isnothing(dir_save) || isdir(dir_save))
        mkdir(dir_save)
    end
    return dir_save
end

function reactionsystem_to_string(parsed_model_SBML::ModelSBMLSystem,
                                  write_to_file::Bool, path_save_model::String,
                                  model_SBML::ModelSBML)::String

    # ReactionSystem
    combinatoric_ratelaws_arg = parsed_model_SBML.int_stoichiometries ? "true" : "false"
    _species_write = parsed_model_SBML.species
    if isempty(model_SBML.rule_variables)
        sps_arg = "sps"
    elseif parsed_model_SBML.no_species == false
        sps_arg = "[sps; vs]"
    else
        sps_arg = "vs"
        _species_write = replace(_species_write, "sps = Catalyst.@species" => "")
    end
    # Parameters might be an empty set
    if parsed_model_SBML.parameters != "\tps = Catalyst.@parameters "
        _rn_write = "\trn = Catalyst.ReactionSystem(_reactions, t, $sps_arg, ps; name=Symbol(\"" *
                    model_SBML.name *
                    "\"), combinatoric_ratelaws=$combinatoric_ratelaws_arg)"
    else
        _rn_write = "\trn = Catalyst.ReactionSystem(_reactions, t, $sps_arg, Any[]; name=Symbol(\"" *
                    model_SBML.name *
                    "\"), combinatoric_ratelaws=$combinatoric_ratelaws_arg)"
    end

    # Create a function returning the ReactionSystem, specie-map, and parameter-map
    _function_write = "function get_reaction_system(foo)\n\n"
    _function_write *= _species_write * "\n"
    if parsed_model_SBML.variables != "\tvs = ModelingToolkit.@variables"
        _function_write *= parsed_model_SBML.variables * "\n"
    end
    if parsed_model_SBML.parameters != "\tps = Catalyst.@parameters "
        _function_write *= parsed_model_SBML.parameters * "\n\n"
    end
    _function_write *= "\tD = Differential(t)\n\n"
    _function_write *= parsed_model_SBML.reactions * "\n\n"
    _function_write *= _rn_write * "\n\n"
    _function_write *= parsed_model_SBML.specie_map * "\n"
    _function_write *= parsed_model_SBML.parameter_map * "\n"
    _function_write *= "\treturn rn, specie_map, parameter_map\nend"

    # In case user request file to be written
    if write_to_file == true
        open(path_save_model, "w") do f
            write(f, _function_write)
        end
    end

    return _function_write
end
