# SBMLImporter.jl 

SBMLImporter.jl is a Julia package for importing dynamic [Systems Biology Markup Language (SBML)](https://sbml.org/) models into either a `JumpProblem` for Gillespie simulations, a `SDEProblem` for chemical Langevin simulations, or an `ODEProblem` for deterministic simulations.

## Major highlights

* Import of SBML models into a [Catalyst.jl](https://github.com/SciML/Catalyst.jl) `ReactionSystem`. This allows for easy conversion to a `JumpProblem`, a `SDEProblem`, or an `ODEProblem`.
* Support for a majority of SBML features, such as dynamic compartments, events, rules, piecewise (`ifelse`) expressions, and units. An extensive feature list can be found [here](@ref support).
* Thoroughly tested against both the [SBML test suite](https://github.com/sbmlteam/sbml-test-suite) and a large collection of published models.
* Integrates with [PEtab.jl](https://github.com/sebapersson/PEtab.jl) for fitting SBML models to data.

!!! note "Star us on GitHub!"
    If you find the package useful in your work please consider giving us a star on [GitHub](https://github.com/sebapersson/SBMLImporter.jl). This will help us secure funding in the future to continue maintaining the package.

## Installation

To install SBMLImporter.jl in the Julia REPL enter

```julia
julia> ] add SBMLImporter
```

or alternatively

```julia
julia> using Pkg; Pkg.add("SBMLImporter")
```

SBMLImporter is compatible with Julia version 1.10 and above. For best performance we strongly recommend using the latest Julia version.

## Getting help

If you have any problems using SBMLImporter, here are some relevant actions:

* Read the [FAQ](@ref FAQ) section in the online documentation.
* Post your questions in the `#sciml-sysbio` channel on the [Julia Slack](https://julialang.org/slack/).
* If you believe you have encountered unexpected behavior or a bug, please open an issue on [GitHub](https://github.com/sebapersson/SBMLImporter.jl).
