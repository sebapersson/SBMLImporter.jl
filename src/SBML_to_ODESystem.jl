"""
    SBML_to_ODESystem(path_SBML::AbstractString;
                      ifelse_to_callback::Bool=true,
                      inline_assignment_rules::Bool=false,
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
- `inline_assignment_rules=true`: Whether to inline assignment rules into model equations. Recomended for model import speed, 
    however, note that it will not be possible to access the rule-variable then via sol[:var]
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
                           inline_assignment_rules::Bool=true,
                           write_to_file::Bool=false,
                           verbose::Bool=true,
                           ret_all::Bool=false,
                           model_as_string::Bool=false) where T <: AbstractString

    rn, specie_map, parameter_map, cb = SBML_to_ReactionSystem(path_SBML; ifelse_to_callback=ifelse_to_callback, 
                                                               inline_assignment_rules=inline_assignment_rules, 
                                                               write_to_file=write_to_file, 
                                                               verbose=verbose, ret_all=true, 
                                                               model_as_string=model_as_string)
    sys = structural_simplify(convert(ODESystem, rn))

    if ret_all == true
        return sys, specie_map, parameter_map, cb
    end

    # Depending on model return what is needed to perform forward simulations
    if isempty(model_SBML.events) && isempty(model_SBML.ifelse_bool_expressions)
        return sys, specie_map, parameter_map
    end
    
    verbose && @info "SBML model with events - output returned as odesys, specie_map, parameter_map, cbset\nFor how to simulate model see documentation"
    return sys, specie_map, parameter_map, cb
end
