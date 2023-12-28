# SBMLImporter.jl
*Julia Importer for Dynamic Models in the SBML Format*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/dev/)
[![Build Status](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml?query=branch%3Amain)

SBMLImporter.jl is a Julia importer for dynamic models specified in the Systems Biology Markup Language (SBML) into a [Catalyst](https://github.com/SciML/Catalyst.jl) `ReactionSystem` - which allows for model simulations with stochastic (Gillespie, SDE) or deterministic (ODE) simulators. It supports many SBML features, such as events, rate-, assignment-, algebraic-rules, dynamic compartment size, and conversion factors. For a detailed list of supported features, see below.

To perform parameter estimation for a SBML model, see [PEtab.jl](https://github.com/sebapersson/PEtab.jl).

A list of differences compared to [SBMLToolkit](https://github.com/SciML/SBMLToolkit.jl) is provided below. For constraint-based modeling, see [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl).

## Installation

SBMLImporter.jl can be installed via the Julia package manager:

```julia
julia> ] add SBMLImporter
```

Alternatively, you can use:

```julia
julia> using Pkg; Pkg.add("SBMLImporter")
```

## Quick Start

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

Alternatively, the model can be imported directly into an `ODESystem` with:

```julia
sys, specie_map, parameter_map = SBML_to_ODESystem(path_SBML)
```

To import more advanced models with events and/or piecewise (ifelse) expressions, see the [documentation](https://sebapersson.github.io/SBMLImporter.jl/stable/).

## Differences compared to SBMLToolkit

The key differences between SBMLToolkit and SBMLImporter are:

* SBMLToolkit works with (and transforms) species to be in amount. SBMLImporter supports species in amount and/or concentration.

* SBMLToolkit has a cleaner interface, as it performs all model processing via Symbolics.jl.

* SBMLImporter has wider event support, including events with directionality. It further processes events without species in the trigger into a `DiscreteCallback`, making simulations more efficient.

* SBMLImporter rewrites SBML piecewise expressions to callbacks if possible instead of using `ifelse`, this improves integration stability and reduces runtime.

* SBMLImporter has more extensive SBML support, passing more tests in the test-suite. It is further the SBML importer for PEtab.jl, which regularly tests against several published models of various sizes.

## Supported SBML Features

SBMLImporter supports many SBML features for SBML models (level 2 or higher). Currently, excluding FBA models, it successfully passes 1257 out of 1785 test cases. The failed test cases cover features currently not supported. If SBMLImporter lacks support for a feature you would like, please file an issue on GitHub. The features not supported are:

* Delay (creating a delay-differential-equations)
* Events with delay
* Events with priority
* Hierarchical models
* Fast reactions
* Parameter or species names corresponding to Julia constants (`pi`, `NaN`, `true`, `false`)
* Certain uncommon math expressions, such as `lt` with three arguments, `implies` etc...

Import might also fail for complicated nested piecewise expressions inside SBML functions.

## Citation

We will soon publish a paper you can cite if you found SBMLImporter.jl helpful in your work.
