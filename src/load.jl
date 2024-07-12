"""
    load_SBML(path; massaction=false, kwargs...)

Import a [SBML](https://sbml.org/) into a `ParsedReactionNetwork` that can simulated
as a `JumpProblem` (Gillespie), a `SDEProblem` (chemical Langevin), or an `ODEProblem`

The keyword `massaction` should be set to `true` if the model reactions follow chemical
mass action [ADD], which often holds for rule based models generated from packages such as
[BioNetGEn](https://bionetgen.org/). **Note** that `massaction=true` is required for
carrying out efficient JumpProblem (Gillespie) simulations, more details [here](add!)

## Keyword arguments
- `ifelse_to_callback=true`: Whether to rewrite `ifelse` (SBML piecewise) expressions to
    [callbacks](ADD). Recommended as this improves simulation runtime and stabillity.
- `inline_assignment_rules=true`: Whether to inline assignment rules into model equations.
    Recommended for lage models, **however**, if `true` is it not possible to access
    assignment rule variables via `sol[:var]`.
- `write_to_file=false`: Whether to write the parsed SBML model to a Julia file in the same directory as the
    SBML file.
- `model_as_string=false`: Whether or not the model (`path_SBML`) is provided as a string.
- `check_massaction=true`: Whether to check if and rewrite a reaction to mass action if possible.
- `massaction=false`: Determines if reactions are treated as mass-action reactions. This option should
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
function load_SBML(path_SBML::AbstractString; massaction::Bool = false, ifelse_to_callback::Bool = true,
                   inline_assignment_rules::Bool = true, write_to_file::Bool = false,
                   model_as_string::Bool = false)
    model_SBML = parse_SBML(path_SBML; massaction = massaction, ifelse_to_callback = ifelse_to_callback, model_as_string = model_as_string, inline_assignment_rules = inline_assignment_rules)
    model_SBML_sys = _to_system_syntax(model_SBML, inline_assignment_rules; massaction = massaction)
    rn, specie_map, parameter_map = _get_reaction_system(model_SBML_sys, model_SBML.name)

    # Build callback functions
    cbset = create_callbacks(rn, model_SBML, model_SBML.name)

    if write_to_file == true && model_as_string == true
        @warn "If the model is provided as a string we do not support writing it to " *
              "file. Hence write_to_file argument is ignored"
    elseif write_to_file == true
        dirsave = _get_dir_save(write_to_file, model_as_string, path_SBML)
        write_reactionsystem(model_SBML_sys, dirsave, model_SBML)
    end
    parsed_rn = ParsedReactionNetwork(rn, specie_map, parameter_map, nothing, nothing)
    return parsed_rn, cbset
end
