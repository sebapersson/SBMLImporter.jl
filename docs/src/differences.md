# Other SBML related Julia packages

There are currently three additional SBML related packages in Julia. [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl) wraps a subset of the [libSBML](https://sbml.org/software/libsbml/) functionality and is used by SBMLImporter and other SBML related packages for parsing SBML models. [COBREXA.jl](https://github.com/COBREXA/COBREXA.jl) is designed for constraint-based metabolic modeling. Constraint-based models, which are often referred to as flux-balance analysis (FBA) models, are not supported by SBMLImporter.

Lastly, [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl) similarly to SBMLImporter imports dynamic models into a `ReactionSystem`. However, we recommend SBMLImporter due to the following key differences:

* SBMLToolkit works with (and transforms) species to be in amount. SBMLImporter supports species units to be either amount or concentration.

* SBMLImporter has wider event support, including events with directionality. It further processes events without species in the trigger into a `DiscreteCallback`, making simulations more efficient.

* SBMLImporter rewrites SBML piecewise expressions to callbacks if possible instead of using `ifelse`, this improves simulation stability and reduces simulation runtime.
  
* When possible, SBMLImporter converts reactions to so-called `MassActionJump`s. This greatly improve performance for most Jump simulations.

* SBMLImporter has more extensive SBML support, passing more tests in the SBML test-suite. It is further the SBML importer for PEtab.jl, which regularly tests against several published models of various sizes.
