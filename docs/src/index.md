# SBMLImporter.jl Documentation

This is the documentation for SBMLImporter.jl, a Julia package for importing ODE models specified in the Systems Biology Markup Language (SBML). This package offers robust SBML support, covering SBML features such as events, rate, assignment, algebraic rules, dynamic compartment size, and conversion factors. For a list of supported features, see [here](@ref support). For a list ofdifferences compared to [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl), see the [README].

To perform parameter estimation for an SBML model, see [PEtab.jl](https://github.com/sebapersson/PEtab.jl) which for importing SBML models, employs SBMLImporter.

## Tutorial

SBMLImporter is a tool for importing SBML models into a [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl)ODESystem, which with other things allows for symbolic model pre-processing. An ODESystem can easily be converted into an ODEProblem and solved using any ODE solver in [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl). If the model includes events, [callbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) are also generated during the import.

!!! info
    The number of arguments returned by `SBML_to_ODESystem` varies depending on whether the model has events. When importing an SBML model, `SBML_to_ODESystem` will inform about the number of returned arguments.

### Importing a Model Without Events and Piecewise

Importing an SBML model without is straightforward. Given the path to a SBML file do:

```julia
using SBMLImporter
sys, specie_map, parameter_map = SBML_to_ODESystem(path_SBML)
```

Here, `sys` is the `ODESystem`, `specie_map` is a mapping for the initial values, and `parameter_map` is a mapping/values for the model parameters. To simulate the model, create an `ODEProblem`:

```julia
using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Solve ODE with Rodas5P solver
sol = solve(prob, Rodas5P())
```

Setting `jac=true` mean that the Jacobian of the ODE is computed symbolically, which is recommended for performance.

Your explanation is detailed, but I made a few adjustments for clarity and consistency:

### Importing a Model with Events

When importing a SBML model events are rewritten to [callbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/). There are two types of callbacks, ContinuousCallbacks or DiscreteCallbacks. The former use root-finding to identify when the event is triggered, while DiscreteCallbacks solves the ODE up to the event time, apply the event, and then proceeds. Since root-finding can be computationally demanding, SBMLImporter rewrite SBML events into DiscreteCallbacks when possible. To keep track of the discrete event times, the importer returns a function for computing said event times given the model parameters:

```julia
sys, specie_map, parameter_map, cb, get_tstops = SBML_to_ODESystem(path_SBML)
```

Here, `cb` represent the model's events, and `get_tstops` is a function to compute the event times given the model parameters. To simulate the model, do:

```julia
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Compute event times
tstops = get_tstops(prob.p, prob.u0)
sol = solve(prob, Rodas5P(), tstops=tstops, callback=callbacks)
```

### Importing a Model with Time-Dependent Piecewise

Piecewise functions in SBML correspond to the Julia `ifelse` function:

```julia
ifelse(cond, value_true, value_false)
```

Here, if `cond==true`, the statement evaluates to `value_true`. While `ifelse` can be directly encoded in the model, this may decrease performance as a discontinuity is added to the model. Therefore, SBMLImporter attempts to rewrite `ifelse` to callbacks (events). Additionally, as `ifelse` can sometimes be active at time `t0`, SBMLImporter provides a function to adjust ifelse statements active at time zero. To import a model with piecewise, do:

```julia
sys, specie_map, parameter_map, cb, get_tstops, ifelse_t0 = SBML_to_ODESystem(path_SBML)
```

Here, `ifelse_t0` is a vector of functions handling piecewise (ifelse) conditions rewritten to events which are active at time zero. To solve the model do:

```julia
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
tstops = get_tstops(prob.p, prob.u0)
# Adjust ifelse statements active at time zero
for _f! in ifelse_t0
    _f!(prob.u0, prob.p)
end
sol = solve(prob, Rodas5P(), tstops=tstops, callback=cb)
```

To prevent the rewriting of `ifelse` to events when creating the ODESystem, set `ifelse_to_callback=false` when calling `SBML_to_ODESystem`.

## Citation

We will soon publish a paper you can cite if you found SBMLImporter.jl helpful in your work.
