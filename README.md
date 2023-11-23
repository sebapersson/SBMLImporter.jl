# SBMLImporter.jl

*Julia Importer for ODE Models in SBML Format*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/dev/)
[![Build Status](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml?query=branch%3Amain)

SBMLImporter.jl is a Julia package forimport Ordinary Differential Equation (ODE) models specified in the Systems Biology Markup Language (SBML) at level 2 or higher. It supports many SBML features, such as events, rate-, assignment-, algebraic-rules, dynamic compartment size, and conversion factors. For a detailed list of supported features, see below.

In comparison to [SBMLToolkit](https://github.com/SciML/SBMLToolkit.jl), SBMLImporter.jl focuses exclusively on ODE models. A comprehensive comparison against SBMLToolkit is provided below. For constraint-based modeling, see [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl).


## Installation

SBMLImporter.jl can be installed using the Julia package manager:

```julia
julia> ] add SBMLImporter
```

Alternatively, you can use:

```julia
julia> using Pkg; Pkg.add("SBMLImporter")
```

## Quick Start

Importing an SBML modelis straightforward. Given the path to a SBML file do:

```julia
using SBMLImporter
sys, specie_map, parameter_map = SBML_to_ODESystem(path_SBML)
```

Here, `sys` is the `ODESystem`, `specie_map` is a mapping for the initial values, and `parameter_map` is a mapping/values for the model parameters. To simulate the model, construct an `ODEProblem` and solve it using an ODE solver:

```julia
using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
sol = solve(prob, Rodas5P())
```

For more advanced models with events and/or piecewise (ifelse) functions, see the documentation.

## Differences with SBMLToolkit

The key differences between SBMLToolkit and SBMLImporter are:

* SBMLToolkit imports models as `ReactionSystems`, allowing simulation with the Gillespie algorithm or transformation into the Langevin SDE. However, it only works with (and transforms) species to be in amount. SBMLImporter transforms models only into `ODESystems`, supporting species in units as amount and concentration.

* SBMLToolkit has a cleaner interface, as it performs all model processing via Symbolics.jl.

* SBMLImporter has broader event support, including events with directionality. It processes events without species in the trigger into `DiscreteCallbacks`, making simulations more efficient.

* SBMLImporter rewrites SBML piecewise expressions to callbacks if possible instead of rewriting to `ifelse`, this improves integration stability and reduces runtime.

* SBMLImporter has more extensive SBML support, passing more tests in the test-suite. It is also the SBML importer for SBMLImporter.jl, which regularly tests against several published SBML models of various sizes.

## Supported SBML Features

SBMLImporter has extensive SBML support and currently passes ADD test cases (excluding hierarchical models and FBA models). If SBMLImporter lacks support for a feature you would like, please file an issue on GitHub. Features currently not supported include:

* Events with priority
* Events with delay
* Delay (creating a delay-differential-equations)
* Fast reactions
* Parameter or species names corresponding to Julia constants (pi, NaN, true, false)

Import might also fail for complicated nested piecewise expressions inside functions. 