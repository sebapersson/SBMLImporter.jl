# SBMLImporter.jl

*Julia Importer for ODE Models in SBML Format*

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

## Quick start

Importing a SBML without events and piecewise functions is straightforward. Given the path to a SBML file (`SBML_to_ODESystem`) do:

```julia
using SBMLImporter
sys, specie_map, parameter_map = SBML_to_ODESystem(path_SBML)
```

Here `sys` is the `ODESystem`, `specie_map` mapping for the initial value, and `parameter_map` mapping/values for the model parameters. Given these components simulating the model is straighforward via building a an `ODEProblem`:

```julia
using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(ode_system, specie_map, tspan, parameter_map, jac=true)
# Solve ODE with Rodas5P solver
sol = solve(prob, Rodas5P())
```

For importing and simulating more advanced models with events and piecewise (ifelse), see the documentation.

## Differences with SBMLToolkit

Here are some key differences between SBMLToolkit and SBMLImporter that can help you choose SBML importer:

* SBMLToolkit imports models as `ReactionSystems`, allowing simulation with the Gillespie algorithm or transformation into the Langevin SDE. However, it only works with (and transforms) species to be in amount. SBMLImporter transforms models only into `ODESystems`, supporting species in units as amount and concentration.

* SBMLImporter has broader event support, including events with directionality. It processes events without species in the trigger into `DiscreteCallbacks`, making simulations more efficient.

* SBMLImporter rewrites SBML piecewise expressions to callbacks if possible instead of rewriting to `ifelse`, this improves integration stability and reduces runtime.

* SBMLImporter has more extensive SBML support, passing more tests in the test-suite. It is also the SBML importer for PEtab.jl, which regularly tests against several published SBML models of various sizes.

## Supported SBML Features

SBMLImporter has extensive SBML support and currently passes ADD test cases (excluding hierarchical models and FBA models). If SBMLImporter lacks support for a feature you would like, please file an issue on GitHub. Features currently not supported include:

* Events with priority
* Events with delay
* Delay (creating a delay-differential-equations)
* Fast reactions
* Parameter or species names corresponding to Julia constants (pi, NaN, true, false)

Import might also fail for complicated nested piecewise expressions inside functions. 