"""
    SBML_to_ODESystem(path_SBML::AbstractString;
                      ifelse_to_callback::Bool=true,
                      write_to_file::Bool=false,
                      verbose::Bool=true,
                      return_all::Bool=false,
                      model_as_string::Bool=false)

Parse an SBML model into a ModelingToolkit `ODESystem` and potentially convert events/piecewise to callbacks.

By default, `structurally_simplified` is called on the `ODESystem` before it is returned.

For information on simulating the `ODESystem`, refer to the documentation.

For testing `path_SBML` can be the model as a string if `model_as_string=true`

!!! note
    The number of returned arguments depends on whether the SBML model has events and/or piecewise expressions (see below).

## Arguments
- `path_SBML`: File path to a valid SBML file (level 2 or higher).
- `ifelse_to_callback=true`: Whether to rewrite `ifelse` (piecewise) expressions to callbacks; recommended for performance.
- `write_to_file=false`: Whether to write the parsed SBML model to a Julia file in the same directory as the SBML file.
- `verbose=true`: Whether or not to display information on the number of return arguments.
- `return_all=true`: Whether or not to return all possible arguments (see below), regardless of whether the model has events.
- `model_as_string=false` : Whether or not the model (`path_SBML`) is provided as str, mainly for testing.

## Returns
- `ode_system`: A ModelingToolkit `ODESystem` that can be converted into an `ODEProblem` and solved.
- `specie_map`: A species map setting initial values; together with the `ODESystem`, it can be converted into an `ODEProblem`.
- `parameter_map` A parameter map setting parameter values; together with the `ODESystem`, it can be converted into an `ODEProblem`.
- `cbset` - **only for models with events/piecewise expressions**: Callbackset (events) for the model.
- `get_tstops`- **Only for models with events/piecewise expressions**: Function computing time stops for discrete callbacks in the `cbset`.

## Examples
```julia
# Import and simulate model without events
using SBMLImporter
sys, specie_map, parameter_map = SBML_to_ODESystem(path_SBML)

using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Solve ODE with Rodas5P solver
sol = solve(prob, Rodas5P())
```
```julia
# Import a model with events
using SBMLImporter
sys, specie_map, parameter_map, cb, get_tstops = SBML_to_ODESystem(path_SBML)

using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Compute event times
tstops = get_tstops(prob.u0, prob.p)
sol = solve(prob, Rodas5P(), tstops=tstops, callback=callbacks)
```
"""
function SBML_to_ODESystem(path_SBML::T;
                           ifelse_to_callback::Bool=true,
                           write_to_file::Bool=false,
                           verbose::Bool=true,
                           ret_all::Bool=false,
                           model_as_string::Bool=false) where T <: AbstractString

    # Intermediate model representation of a SBML model which can be processed into
    # an ODESystem
    model_SBML = build_SBML_model(path_SBML; ifelse_to_callback=ifelse_to_callback, model_as_string=model_as_string)
    rule_variables = unique(vcat(model_SBML.rate_rule_variables, model_SBML.assignment_rule_variables,
                                 model_SBML.algebraic_rule_variables))

    # If model is written to file save it in the same directory as the SBML-file
    if model_as_string == false
        dir_save = joinpath(splitdir(path_SBML)[1], "SBML")
    else
        dir_save = joinpath(@__DIR__, "SBML")
    end
    if write_to_file == true && !isdir(dir_save)
        mkdir(dir_save)
    end
    path_save_model = joinpath(dir_save, model_SBML.name * ".jl")

    # Build the ReactionSystem. Must be done via Meta-parse, because if a function is used 
    # via @RuntimeGeneratedFunction runtime is very slow for large models
    (_species_write, _specie_map_write, _variables_write, 
     _parameters_write, _parameter_map_write, _reactions_write, 
     no_species, integer_stoichiometries) = _reactionsystem_from_SBML(model_SBML, rule_variables)
    # The model can have have only species or only variables or both. If it has variables they 
    # are given via SBML rules
    eval(Meta.parse("ModelingToolkit.@variables t"))
    eval(Meta.parse("Differential(t)"))
    sps = eval(Meta.parse(split(_species_write, "\n")[2]))
    vs = isempty(rule_variables) ? Any[] : eval(Meta.parse(_variables_write))
    if isempty(rule_variables)
        sps_arg = sps
    elseif no_species == false
        sps_arg = [sps; vs]
    else
        sps_arg = vs
    end                                  
    # Parameters can not be an empty collection
    if _parameters_write != "\tps = Catalyst.@parameters "
        ps = eval(Meta.parse(_parameters_write))
    else
        ps = Any[]
    end
    _reactions = eval(Meta.parse(_reactions_write))
    combinatoric_ratelaws = integer_stoichiometries ? true : false
    # Build reaction system from its components
    reaction_system = Catalyst.ReactionSystem(_reactions, t, sps_arg, ps; name=Symbol(model_SBML.name), combinatoric_ratelaws=combinatoric_ratelaws)
    specie_map = eval(Meta.parse(_specie_map_write))
    parameter_map = eval(Meta.parse(_parameter_map_write))
    
    # Build callback functions
    cbset, callback_str = create_callbacks_SBML(reaction_system, model_SBML, model_SBML.name)
    # Convert into an ODESystem
    if isempty(model_SBML.algebraic_rules)
        ode_system = structural_simplify(convert(ODESystem, reaction_system))
    else
        ode_system = structural_simplify(dae_index_lowering(convert(ODESystem, reaction_system)))
    end

    # if model is written to file write the callback
    if write_to_file == true
        path_save = joinpath(dir_save, model_SBML.name * "_callbacks.jl")
        io = open(path_save, "w")
        write(io, callback_str)
        close(io)
        _ = reactionsystem_to_string(_species_write, _specie_map_write, _variables_write, 
                                     _parameters_write, _parameter_map_write, _reactions_write, 
                                     no_species, integer_stoichiometries, write_to_file, path_save_model, 
                                     rule_variables, model_SBML)
    end

    if ret_all == true
        return ode_system, specie_map, parameter_map, cbset
    end

    # Depending on model return what is needed to perform forward simulations
    if isempty(model_SBML.events) && isempty(model_SBML.ifelse_bool_expressions)
        return ode_system, specie_map, parameter_map
    end

    verbose && @info "SBML model with events - output returned as odesys, specie_map, parameter_map, cbset\nFor how to simulate model see documentation"
    return ode_system, specie_map, parameter_map, cbset
end
