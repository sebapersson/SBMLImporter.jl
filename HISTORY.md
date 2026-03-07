# Breaking updates and feature summaries across releases

## SBMLImporter 4.0

SBMLImporter 4.0 is a breaking release triggered by the Catalyst.jl v16 update. Major
changes are:

- Updated Catalyst to v16. Since Catalyst now uses a more flexible metadata interface,
  `load_SBML` now returns a `ReactionSystem` (rather than a `ParsedReactionSystem`). The
  initial value map and parameter map are stored as metadata on the returned system. See
  below for how to import a model with `load_SBML`.
- Dropped ModelingToolkit as a direct dependency and switched to ModelingToolkitBase. Since
  ModelingToolkit v11+ is AGPL-licensed, ModelingToolkit is no longer a direct dependency in
  order to keep SBMLImporter.jl and its dependencies permissively licensed.

### Updated `load_SBML` return type and usage

Models are now imported as a `ReactionSystem`, where the initial value map (`u0_map`) and
parameter value map (`parameter_map`) are stored as metadata. Importing as an `ODEProblem`
now looks like:

```julia
# New syntax
using SBMLImporter, ModelingToolkitBase, OrdinaryDiffEq
rn, cb = load_SBML(path_sbml)
u0 = get_u0_map(rn)
ps = get_parameter_map(rn)
sys = mtkcompile(ode_model(rn))
tspan = (0.0, 10.0)
oprob = ODEProblem(sys, merge(Dict(u0), Dict(ps)), tspan; jac = true)

# Old syntax
using SBMLImporter, ModelingToolkit, OrdinaryDiffEq
prn, cb = load_SBML(path_sbml)
sys = structural_simplify(convert(ODESystem, prn.rn))
tspan = (0.0, 10.0)
oprob = ODEProblem(sys, prn.u0, tspan, prn.p; jac = true)
```

Examples for importing into a `JumpProblem` or an `SDEProblem` are available in the online
documentation.

### ModelingToolkitBase and licensing

Since ModelingToolkit v11+ is AGPL-licensed (see discussion
[here](https://discourse.julialang.org/t/modelingtoolkit-v11-library-split-and-licensing-community-feedback-requested/134396)),
SBMLImporter.jl now depends on the MIT-licensed ModelingToolkitBase instead. SBMLImporter.jl
aims to remain permissively licensed, and with this update none of its direct dependencies
are AGPL-licensed (or more restrictive).

Practically, this has limited impact. The imported `ReactionSystem` can still be used with
ModelingToolkit for downstream tasks (e.g. rewriting SBML models with algebraic rules).

## SBMLImporter 3.2

Added a flag to prevent SBML events with `initialValue=false` from triggering at `t0`. This
is needed to properly support PEtab v2 time courses.

## SBMLImporter 3.1

Added support for Julia 1.12.

## SBMLImporter 3.0

SBMLImporter 3.0 is a breaking release triggered by Catalyst.jl updating to v15. The
SBMLImporter API is unchanged, but as SBML models are imported into a
`Catalyst.ReactionSystem`, only Catalyst v15 features are supported for downstream modeling
tasks promptig the major update.

## SBMLImporter 2.8

Added support for 5-argument `piecewise` functions.

## SBMLImporter 2.7

- Added support for importing SBML models as an `ODEProblem` (required for PEtab SciML
  integration).
- Bumped minimum supported Julia version to 1.10.

## SBMLImporter 2.6

Updated internal code to use the new `t` and `time_deriv` functions with ModelingToolkit v9
(no user-facing API changes expected).

## SBMLImporter 2.5

Updated internal event-time handling to use new Symbolics solvers (no user-facing API
changes expected).

## SBMLImporter 2.4

Fixed internal issues affecting PEtab.jl integration (no user-facing API changes expected).

## SBMLImporter 2.3

Fixed internal issues affecting PEtab.jl integration (no user-facing API changes expected).

## SBMLImporter 2.2

Added compartment metadata. For each model species, its associated compartment can be
retrieved via:

```julia
using SBMLImporter, Catalyst
prn, cb = load_SBML(path_SBML)
sbml_species = species(prn.rn)
c = getcompartment(sbml_species[1])
```

## SBMLImporter 2.1

- Added SBML reaction id and name as metadata for Catalyst reactions.
- Added an option to `load_SBML` to avoid inlining kinetic-law parameters.

## SBMLImporter 2.0

Breaking release due to updating to ModelingToolkit v9. Versions prior to 2.0 are not
supported.
