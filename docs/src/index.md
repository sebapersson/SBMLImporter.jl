# SBMLImporter.jl 

SBMLImporter.jl is an importer for dynamic models defined in the Systems Biology Markup Language (SBML). It supports most SBML features, such as events, dynamic compartment sizes, and rate, assignment, and algebraic rules. For a complete list of supported features see [here](@ref support). For differences compared to [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl) see the [README](https://github.com/sebapersson/SBMLImporter.jl).

To perform parameter estimation for a SBML model, see [PEtab.jl](https://github.com/sebapersson/PEtab.jl).

## Installation

To install SBMLImporter.jl in the Julia REPL enter

```julia
julia> ] add SBMLImporter
```

or alternatively

```julia
julia> using Pkg; Pkg.add("SBMLImporter")
```

SBMLImporter.jl is compatible with Julia version 1.6 and above. For best performance we strongly recommend using Julia version 1.10.

## Tutorial

SBMLImporter import SBML models as a [Catalyst](https://github.com/SciML/Catalyst.jl) `ReactionSystem`. This provides several benefits, such as symbolic model pre-processing for efficient simulations. The imported `ReactionSystem` can be converted to a `JumpProblem` for Gillespie simulations, a `SDEProblem` for Langevin SDE simulations, or an `ODEProblem` for deterministic ODE simulations

As example, consider the Brusselator model (the SBML file can be downloaded from [here](https://github.com/sebapersson/SBMLImporter.jl/blob/main/test/Models/brusselator.xml)). The first step is to import the model with `load_SBML`:

```@example 1
using SBMLImporter
path_SBML = joinpath(@__DIR__, "..", "..", "test", "Models", "brusselator.xml") # hide
prnbng, cb = load_SBML(path_SBML)
nothing # hide
```

This returns two outputs a `ParsedReactionSystem` (`prnbng`) and a `CallbackSet` (`cb`). The `ParsedReactionSystem` includes the reaction system (`prnbng.rn`), a map for the initial values of each species (`prnbng.u₀`), and a map setting the model parameter values (`prnbng.p`). The `CallbackSet` holds any potential SBML events, along with SBML piecewise functions that have been parsed into events.

### Gillespie simulations

To perform Gillespie simulations, convert the reaction-system `prnbng.rn` into a `JumpProblem`.

```@example 1
using JumpProcesses
using Random # hide
Random.seed!(1) # hide
tspan = (0.0, 10.0)
dprob = DiscreteProblem(prnbng.rn, prnbng.u₀, tspan, prnbng.p)
jprob = JumpProblem(prnbng.rn, dprob, Direct())
nothing # hide
```

The `JumpProblem` can be solved with any solver from the [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) package, such as the `SSAStepper`:

```@example 1
using Plots
sol = solve(jprob, SSAStepper(), callback=cb)
plot(sol; lw=2)
```

For more information on Gillespie simulations, see the documentation for [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl).

!!! warn
    For efficient Gillespie simulations two conditions must be met: the model should be a mass-action model and each species should have units amount. This translates to ensuring that every species has the attribute `hasOnlySubstanceUnits=true`, and no rule variables are used in the kinetic math expressions for the SBML reactions.

### SDE simulations

To perform SDE simulations, convert the reaction-system `prnbng.rn` into a `SDEProblem`.

```@example 1
using StochasticDiffEq
tspan = (0.0, 10.0)
sprob = SDEProblem(prnbng.rn, prnbng.u₀, tspan, prnbng.p)
nothing # hide
```

The `SDEProblem` can be solved with any solver from the [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) package, such as the `LambaEM` solver:

```@example 1
sol = solve(sprob, LambaEM(), callback=cb)
plot(sol; lw=2)
```

For more information on SDE simulations, see the documentation for [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

### ODE simulations

To perform ODE simulations, convert the reaction-system `prnbng.rn` into an `ODEProblem`.

```@example 1
using ModelingToolkit, OrdinaryDiffEq
tspan = (0.0, 10.0)
sys = convert(ODESystem, prnbng.rn)
oprob = ODEProblem(sys, prnbng.u₀, tspan, prnbng.p, jac=true)
nothing # hide
```

Here `jac=true` means that the ODE Jacobian is computed symbolically which can help with simulation performance. The `ODEProblem` can be solved with any solver from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package, such as the `Rodas5` solver:

```@example 1
sol = solve(oprob, Rodas5(), callback=cb)
plot(sol; lw=2)
```

For more information on ODE simulations, see the documentation for [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

## Citation

We will soon publish a paper you can cite if you found SBMLImporter.jl helpful in your work.
