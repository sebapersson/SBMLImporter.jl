# SBMLImporter.jl

This is the documentation for SBMLImporter.jl, a Julia importer for dynamic models specified in the Systems Biology Markup Language (SBML). This importer supports many SBML features such as events, dynamic compartments size, rate-, assignment-, and algebraic-rules. For a list of supported features, see [here](@ref support). For a list of differences compared to [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl), see the [README](https://github.com/sebapersson/SBMLImporter.jl).

To perform parameter estimation for a SBML model, see [PEtab.jl](https://github.com/sebapersson/PEtab.jl).

## Tutorial

SBMLImporter is a tool for importing SBML models into a [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) `ODESystem` or a [Catalyst](https://github.com/SciML/Catalyst.jl) `ReactionSystem`. This offers several benefits, such as symbolic model pre-processing for efficient simulations. An `ODESystem` can easily be converted into an `ODEProblem` and solved using any ODE solver in [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl), while a `ReactionSystem` can for example easily be converted into an `ODESystem` or `SDESystem`. If the model includes events, [callbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) are generated during the import.

!!! note
    The number of arguments returned by `SBML_to_ReactionSystem` and `SBML_to_ODESystem` varies depending on whether the model has events. When importing an SBML model, the import function will inform about the number of returned arguments.

### Importing a Model Without Events and Without Piecewise Expressions

Importing an SBML model is straightforward. Given the path to a SBML file to import into a `ReactionSystem` do:

```julia
using SBMLImporter
rn, specie_map, parameter_map = SBML_to_ReactionSystem(path_SBML)
```

Here, `rn` is the `ReactionSystem` that for example can be converted into an `ODESystem` or a `SDESystem`, `specie_map` is a mapping for the initial values, and `parameter_map` is a mapping/values for the model parameters. To simulate the model with an ODE-solver, construct an `ODEProblem` and solve it using any ODE solver from [OrdinaryDiffeq](https://github.com/SciML/OrdinaryDiffEq.jl):

```julia
using OrdinaryDiffEq
sys = convert(ODESystem, rn)
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
sol = solve(prob, Rodas5P())
```

Setting `jac=true` mean that the Jacobian of the ODE is computed symbolically, which is recommended for performance. To get the order of the species and parameters in the model do:

```julia
using ModelingToolkit
states(sys) # species
parameters(sys)
```

Alternatively, the model can be imported directly into an `ODESystem` with:

```julia
sys, specie_map, parameter_map = SBML_to_ODESystem(path_SBML)
```

From this point the documentation focuses on ODE-models, but any model that can be imported as an `ODESystem` can also be imported as a `ReactionSystem`.

### Importing a Model with Events

When importing a SBML model with events, the events are rewritten to [callbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/). There are two types of callbacks, `ContinuousCallback` and `DiscreteCallback`. The former use root-finding to identify when the event is triggered, while a `DiscreteCallback` solves the ODE up to the event time, applies the event, and then proceeds. Since root-finding can be computationally demanding, SBMLImporter rewrites a SBML events into a `DiscreteCallback` when possible. To keep track of the discrete event times, the importer also returns a function for computing event times given the model parameters:

```julia
sys, specie_map, parameter_map, cb, get_tstops = SBML_to_ODESystem(path_SBML)
```

Here, `cb` represent the model's events, and `get_tstops` is a function to compute the event times. To simulate the model, do:

```julia
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Compute event times
tstops = get_tstops(prob.u0, prob.p)
sol = solve(prob, Rodas5P(), tstops=tstops, callback=callbacks)
```

### Importing a Model with Time-Dependent Piecewise Expressions

In SBML Piecewise expressions correspond to the Julia `ifelse` function:

```julia
ifelse(cond, value_true, value_false)
```

If `cond==true`, the statement evaluates to `value_true`. While `ifelse` statements can be directly encoded in the model, this may decrease performance as a discontinuity is added. Therefore, SBMLImporter attempts to rewrite `ifelse` to callbacks (events). Additionally, as `ifelse` can sometimes be active at time `t0`, SBMLImporter provides a function to adjust ifelse rewritten callbacks at time zero:

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

To not rewrite `ifelse` to events when creating the ODESystem, set `ifelse_to_callback=false` when calling `SBML_to_ODESystem`.

## Citation

We will soon publish a paper you can cite if you found SBMLImporter.jl helpful in your work.
