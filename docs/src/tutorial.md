# Tutorial

This overarching tutorial of SBMLImporter will cover how to import an SBML model into a `JumpProblem` for Gillespie simulations, a `SDEProblem` for chemical Langevin simulations, or an `ODEProblem` for deterministic simulations.

## Input - a valid SBML model file

SBMLImporter only requires one input: a valid SBML file. SBML files can be created from various sources, such as the graphical interface of [COPASI](https://copasi.org/), by converting rule-based [BioNetGen](https://github.com/RuleWorld/bionetgen) models, or by using any other [SBML exporting tools](https://sbml.org/software/). Additionally, a large collection of published SBML models are hosted on [BioModels](https://www.ebi.ac.uk/biomodels/).

In this tutorial, we will use the [Brusselator model](https://en.wikipedia.org/wiki/Brusselator), whose SBML file can be downloaded from [here](https://github.com/sebapersson/SBMLImporter.jl/blob/main/test/Models/brusselator.xml).

## Importing a model

To import an SBML model use `load_SBML`:

```@example 1
using SBMLImporter
path = joinpath(@__DIR__, "..", "..", "test", "Models", "brusselator.xml") # hide
prn, cb = load_SBML(path; massaction=true)
nothing # hide
```

Here, `massaction=true` informs the importer that the model follows chemical [mass action](https://en.wikipedia.org/wiki/Law_of_mass_action) kinetics. This is required for efficient jump (Gillespie type) simulations.

`load_SBML` returns two outputs: a `ParsedReactionSystem` (`prn`) and a `CallbackSet` (`cb`). The `ParsedReactionSystem` includes the reaction system (`prn.rn`), a map for the initial condition values of each species (`prn.u0`), and a map for setting the model parameter values (`prn.p`). The `CallbackSet` holds any potential SBML events, along with SBML piecewise functions that have been parsed into events.

!!! info
    The `massaction` argument is only relevant for jump simulations. If you import the model as a `SDEProblem` or an `ODEProblem`, you can (and should) ignore this keyword argument.

## Jump simulations

To perform jump simulations (e.g. using Gillespie's algorithm), convert the reaction system (`prn.rn`) into a `JumpProblem` (which requires first creating an intermediary `DiscreteProblem`).

```@example 1
using JumpProcesses
using Random # hide
Random.seed!(1) # hide
tspan = (0.0, 10.0)
dprob = DiscreteProblem(prn.rn, prn.u0, tspan, prn.p)
jprob = JumpProblem(prn.rn, dprob, Direct())
nothing # hide
```

The `JumpProblem` can be simulated with any solver from the [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) package, such as the `SSAStepper`:

```@example 1
using Plots
sol = solve(jprob, SSAStepper(), callback=cb)
default(left_margin=12.5Plots.Measures.mm, bottom_margin=12.5Plots.Measures.mm) # hide
plot_args = (titlefontsize=12, guidefontsize=12, tickfontsize=12, legendfontsize=12, size=(600*1.2, 400*1.2), lw = 3, xlabel = "Time [s]", ylabel = "Protein amount") # hide
plot(sol; lw = 3, xlabel = "Time [s]", ylabel = "Protein amount")
plot(sol; plot_args...) # hide
```

For more information on jump simulations, see [JumpProcesses.jl's documentation](https://github.com/SciML/JumpProcesses.jl).

## SDE simulations

To perform chemical Langevin SDE simulations, convert the reaction-system `prn.rn` into a `SDEProblem`.

```@example 1
using StochasticDiffEq
tspan = (0.0, 10.0)
sprob = SDEProblem(prn.rn, prn.u0, tspan, prn.p)
nothing # hide
```

The `SDEProblem` can be simulated with any solver from the [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) package, such as the `LambaEM` solver:

```@example 1
sol = solve(sprob, LambaEM(), callback=cb)
plot(sol; lw = 3, xlabel = "Time [s]", ylabel = "Protein amount")
plot(sol; plot_args...) # hide
```

For more information on SDE simulations, see [StochasticDiffEq.jl's documentation](https://github.com/SciML/StochasticDiffEq.jl).

## ODE simulations

To perform ODE simulations, convert the reaction system (`prn.rn`) into an `ODEProblem`.

```@example 1
using ModelingToolkit, OrdinaryDiffEq
tspan = (0.0, 10.0)
sys = convert(ODESystem, prn.rn)
oprob = ODEProblem(sys, prn.u0, tspan, prn.p, jac=true)
nothing # hide
```

Here `jac=true` means that the ODE Jacobian is computed symbolically which can help with simulation performance. The `ODEProblem` can be simulated with any solver from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package, such as the `Rodas5P` solver:

```@example 1
sol = solve(oprob, Rodas5P(), callback=cb)
plot(sol; lw = 3, xlabel = "Time [s]", ylabel = "Protein amount")
plot(sol; plot_args...) # hide
```

For more information on ODE simulations, see [OrdinaryDiffEq.jl's documentation](https://github.com/SciML/OrdinaryDiffEq.jl).

## Next Steps

For additional modeling tasks that can be carried out with a Catalyst `ReactionSystem` (e.g., bifurcation analysis, parameter estimation, etc.), see the Catalyst [documentation](https://github.com/SciML/Catalyst.jl). If you encounter any problems with importing an SBML model, first consult the [FAQ](@ref FAQ). Additional details on importer options can be found in the [API](@ref API).
