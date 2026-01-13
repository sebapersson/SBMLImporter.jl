# Tutorial

This tutorial demonstrates how to import an SBML model with SBMLImporter and simulate it as
a `JumpProblem` (Gillespie/SSA), `SDEProblem` (chemical Langevin), or `ODEProblem`
(deterministic ODEs). It also describes how to modify parameters and initial conditions, as
we as how to index solution trajectories by SBML IDs.

## Input: a valid SBML file

SBMLImporter only requires one input: a valid SBML file. SBML files can be created from
various sources, such as the graphical interface of [COPASI](https://copasi.org/), by
converting rule-based [BioNetGen](https://github.com/RuleWorld/bionetgen) models, or by
using any other [SBML exporting tools](https://sbml.org/software/). Additionally, a large
collection of published SBML models are hosted on
[BioModels](https://www.ebi.ac.uk/biomodels/).

This tutorial will use the
[Brusselator model](https://en.wikipedia.org/wiki/Brusselator), whose SBML file can be
downloaded from this
[link](https://github.com/sebapersson/SBMLImporter.jl/blob/main/test/Models/brusselator.xml).

## Importing a model

SBML models are imported with `load_SBML`:

```@example 1
using SBMLImporter
path_sbml = joinpath(@__DIR__, "..", "..", "test", "Models", "brusselator.xml") # hide
prn, cb = load_SBML(path_sbml; massaction=true)
nothing # hide
```

The keyword `massaction=true` indicates that reactions follow chemical
[mass-action kinetics](https://en.wikipedia.org/wiki/Law_of_mass_action). This enables
efficient jump (SSA-type) simulation.

`load_SBML` returns two outputs: a `ParsedReactionSystem` (`prn`) and a `CallbackSet`
(`cb`). The `CallbackSet` holds any SBML events, as well as any `piecewise` functions
parsed into events. The `ParsedReactionSystem` includes a [Catalyst](https://github.com
SciML/Catalyst.jl) `ReactionSystem` (`prn.rn`), a map for the initial values of each
species (`prn.u0`), and a map for setting the model parameter values (`prn.p`). Many
modeling tasks can be performed with a Catalyst `ReactionSystem`. For example,
model reactions can be inspected with:

```@example 1
using Catalyst
reactions(prn.rn)
prn.rn # hide
```

::: info

The `massaction` keyword only affects jump simulations. For `SDEProblem` and
`ODEProblem` workflows, it can be ignored.

:::

## Jump simulations

Jump simulations (e.g. Gillespie/SSA) are run by constructing a `JumpProblem`. This
requires an intermediate `DiscreteProblem`:

```@example 1
using JumpProcesses
using Random # hide
Random.seed!(1) # hide
tspan = (0.0, 10.0)
dprob = DiscreteProblem(prn.rn, prn.u0, tspan, prn.p)
jprob = JumpProblem(prn.rn, dprob, Direct())
nothing # hide
```

The resulting `JumpProblem` can be simulated with solvers from
[JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl), for example `SSAStepper`:

```@example 1
using Plots
sol = solve(jprob, SSAStepper(), callback=cb)
default(left_margin=12.5Plots.Measures.mm, bottom_margin=12.5Plots.Measures.mm, size = (600*1.1, 400 * 1.1), palette = ["#CC79A7", "#009E73", "#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#F0E442"], linewidth=3.0) # hide
plot(sol; xlabel="Time [s]", ylabel="Species amount")
```

Further details are available in the
[JumpProcesses.jl documentation](https://github.com/SciML/JumpProcesses.jl). A guide to
choosing a jump simulation algorithm can be found
[here](https://docs.sciml.ai/JumpProcesses/dev/jump_types#Jump-Aggregators-for-Exact-Simulation).

## SDE simulations

Chemical Langevin simulations are run by constructing an `SDEProblem` from the reaction system:

```@example 1
using StochasticDiffEq
tspan = (0.0, 10.0)
sprob = SDEProblem(prn.rn, prn.u0, tspan, prn.p)
nothing # hide
```

The `SDEProblem` can be simulated with solvers from
[StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl), for example `LambaEM`:

```@example 1
sol = solve(sprob, LambaEM(), callback=cb)
plot(sol; xlabel="Time [s]", ylabel="Species amount")
```

For more information, see the
[StochasticDiffEq.jl documentation](https://github.com/SciML/StochasticDiffEq.jl).

## ODE simulations

Deterministic simulations can be obtained by converting the reaction system into an
`ODESystem`:

```@example 1
using ModelingToolkit
sys = structural_simplify(convert(ODESystem, prn.rn))
nothing # hide
```

Calling `structural_simplify` inlines assignment rules into the final equations. As a
sanity check, the generated ODEs can be inspected:

```@example 1
equations(sys)
```

An `ODEProblem` can then be constructed:

```@example 1
using OrdinaryDiffEq
tspan = (0.0, 10.0)
oprob = ODEProblem(sys, prn.u0, tspan, prn.p, jac=true)
nothing # hide
```

With `jac=true`, the Jacobian is generated symbolically, which often improves performance.
The `ODEProblem` can be simulated with solvers from
[OrdinaryDiffEq.jl](https://github.comSciML/OrdinaryDiffEq.jl), for example `Rodas5P`:

```@example 1
sol = solve(oprob, Rodas5P(), callback=cb)
plot(sol; xlabel="Time [s]", ylabel="Species amount")
```

For more information, see the
[OrdinaryDiffEq.jl documentation](https://github.com/SciML OrdinaryDiffEq.jl). An ODE
solver selection guide is available
[here](https://sebapersson.github.io/PEtab.jl/stable/default_options/).

!!! tip "Importing large models"
    For large models (>1000 species), symbolic Jacobian construction can be time-consuming.
    If import time dominates, setting `jac=false` can help.

## Indexing parameters, species, and solutions

Problems imported by `load_SBML` (`ODEProblem`, `SDEProblem`, or `JumpProblem`) support
[SymbolicIndexingInterface](https://github.com/SciML/SymbolicIndexingInterface.jl),
enabling model variables to be read and updated by SBML IDs.

Parameter values can be modified via `prob.ps`. For example, parameter `A` can be set with:

```@example 1
oprob.ps[:A] = 2.0
nothing # hide
```

Initial conditions can be modified by indexing the problem. For example, the initial value of `X` can be set with:

```@example 1
oprob[:X] = 2.0
nothing # hide
```

After solving, the same symbolic indexing can be used to retrieve trajectories, e.g. for
species `X`:

```@example 1
sol = solve(oprob, Rodas5P(), callback=cb)
sol[:X]
```

The same solution indexing applies to other SBML variables (e.g. assignment rules). If
`inline_assignment_rules=true` was passed to `load_SBML`, assignment rules are inlined into
the model equations during import and cannot be accessed this way.

## Next steps

For downstream modeling tasks supported by Catalyst `ReactionSystem`s (e.g. bifurcation
analysis), see the [Catalyst documentation](https://github.com/SciML/Catalyst.jl). For
parameter estimation of SBML models, see the
[PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/)

If importing an SBML model fails, please consult the [FAQ](@ref FAQ) first. Additional
import options and keyword arguments are documented in the [API](@ref API).
