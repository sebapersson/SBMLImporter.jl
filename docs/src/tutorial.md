# Tutorial

This overarching tutorial of SBMLImporter covers how to import an SBML model into a `JumpProblem` for Gillespie simulations, a `SDEProblem` for chemical Langevin simulations, or an `ODEProblem` for deterministic simulations. It also covers how to change parameter and/or initial values for an imported model.

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

Here, `massaction=true` informs the importer that the model follows chemical [mass action](https://en.wikipedia.org/wiki/Law_of_mass_action) kinetics. This is required for efficient jump (Gillespie type) simulations. Additional keyword arguments can be found in the [API](@ref API) documentation.

`load_SBML` returns two outputs: a `ParsedReactionSystem` (`prn`) and a `CallbackSet` (`cb`). The `CallbackSet` holds any SBML events, as well as any `piecewise` functions parsed into events. The `ParsedReactionSystem` includes a [Catalyst](https://github.com/SciML/Catalyst.jl) `ReactionSystem` (`prn.rn`), a map for the initial values of each species (`prn.u0`), and a map for setting the model parameter values (`prn.p`). Many modeling tasks can be performed with a Catalyst `ReactionSystem`, for example, the system reactions can be inspected by:

```@example 1
using Catalyst
reactions(prn.rn)
prn.rn # hide
```

!!! note
    The `massaction` keyword argument is only relevant for jump simulations. If you import the model as a `SDEProblem` or an `ODEProblem`, you can (and should) ignore this keyword argument.

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
default(left_margin=12.5Plots.Measures.mm, bottom_margin=12.5Plots.Measures.mm, size = (600*1.1, 400 * 1.1), palette = ["#CC79A7", "#009E73", "#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#F0E442"], linewidth=3.0) # hide
plot(sol; xlabel = "Time [s]", ylabel = "Protein amount")
```

For more information on jump simulations, see [JumpProcesses.jl's documentation](https://github.com/SciML/JumpProcesses.jl). A guide to choosing simulation algorithm can be found [here](https://docs.sciml.ai/JumpProcesses/dev/jump_types/#Jump-Aggregators-for-Exact-Simulation).

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
plot(sol; xlabel = "Time [s]", ylabel = "Protein amount")
```

For more information on SDE simulations, see [StochasticDiffEq.jl's documentation](https://github.com/SciML/StochasticDiffEq.jl).

## ODE simulations

To perform ODE simulations, first convert the reaction system (`prn.rn`) into an `ODESystem`.

```@example 1
using ModelingToolkit
sys = structural_simplify(convert(ODESystem, prn.rn))
nothing # hide
```

Here we call `structural_simplify` to inline any potential assignment rules into the model equations. Given an `ODESystem`, a good sanity check is to inspect the generated ODEs:

```@example 1
equations(sys)
```

Once we have an `ODESystem`, we can build an `ODEProblem`:

```@example 1
using OrdinaryDiffEq
tspan = (0.0, 10.0)
oprob = ODEProblem(sys, prn.u0, tspan, prn.p, jac=true)
nothing # hide
```

Setting `jac=true` means that the ODE Jacobian is computed symbolically, which often improves simulation performance. The `ODEProblem` can be simulated with any solver from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package, such as the `Rodas5P` solver:

```@example 1
sol = solve(oprob, Rodas5P(), callback=cb)
plot(sol; xlabel = "Time [s]", ylabel = "Protein amount")
```

For more information on ODE simulations, see [OrdinaryDiffEq.jl's documentation](https://github.com/SciML/OrdinaryDiffEq.jl). A guide to choosing ODE solver can be found [here](https://sebapersson.github.io/PEtab.jl/stable/default_options/).

!!! tip "Importing large models"
    For large models with more than 1000 species, constructing the Jacobian for the `ODEProblem` can be time-consuming. Therefore, if model load time is a problem, setting `jac = false` can help.

## Changing Parameter Values and/or Initial Values

One of the most common modeling operations is to modify model parameter values and/or initial values. Because `SBMLImporter` imports models as a Catalyst `ReactionSystem` that supports the [SymbolicIndexingInterface](https://github.com/SciML/SymbolicIndexingInterface.jl), these modifications are straightforward to perform at the `ODEProblem`, `SDEProblem`, or `JumpProblem` level.

Parameter values can be modified using `prob.ps[:param_id]`. For example, for the Brusselator imported as an `ODEProblem` above, the value of parameter `A` can be changed with:

```@example 1
oprob.ps[:A] = 2.0
nothing # hide
```

Initial values can be modified with `prob[:specie_id]`. For example, to change the initial value of the species `X` to `2.0`, do:

```@example 1
oprob[:X] = 2.0
nothing # hide
```

## Next Steps

For additional modeling tasks that can be carried out with a Catalyst `ReactionSystem` (e.g., bifurcation analysis, parameter estimation, etc.), see the Catalyst [documentation](https://github.com/SciML/Catalyst.jl). If you encounter any problems with importing an SBML model, first consult the [FAQ](@ref FAQ). Additional details on importer options can be found in the [API](@ref API).
