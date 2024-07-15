# SBMLImporter.jl
*Julia Importer for Dynamic Models in the SBML Format*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/dev/)
[![Build Status](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![codecov](https://codecov.io/gh/sebapersson/SBMLImporter.jl/graph/badge.svg?token=1J7RKUPZOM)](https://codecov.io/gh/sebapersson/SBMLImporter.jl)

SBMLImporter.jl is a Julia package for importing dynamic [Systems Biology Markup Language (SBML)](https://sbml.org/) models into either a `JumpProblem` for Gillespie simulations, a `SDEProblem` for chemical Langevin simulations, or an `ODEProblem` for deterministic simulations. Some major highlights of SBMLImporter are:

* It imports models into a [Catalyst.jl](https://github.com/SciML/Catalyst.jl) `ReactionSystem`. This allows for easy conversion to a `JumpProblem`, `SDEProblem`, or `ODEProblem`.
* It supports a majority of SBML features, such as dynamic compartments, events, rules, piecewise (ifelse) expressions, and units.
* It is thoroughly tested against both the SBML test suite and a large collection of published models.
* It integrates with [PEtab.jl](https://github.com/sebapersson/PEtab.jl) for fitting SBML models to data.

Additional information on features and tutorials can be found in the [documentation](https://sebapersson.github.io/SBMLImporter.jl/stable/).

## Related packages

There are currently three SBML related packages in Julia:

1. [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl) wraps a subset of the [libSBML](https://sbml.org/software/libsbml/) functionality and is used by SBMLImporter and other packages for parsing SBML models.
2. [COBREXA.jl](https://github.com/COBREXA/COBREXA.jl) is designed for constraint-based metabolic modeling. Constraint-based, or flux-balance analysis (FBA), is not supported by SBMLImporter.
3. [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl) similarly to SBMLImporter imports dynamic models into a `ReactionSystem`. However, we recommend SBMLImporter as it supports more SBML features, has more efficient event handling, and supports efficient `JumpProblem` (Gillespie) simulations. An extensive list of differences can be found in the documentation.

## Citation

We will soon publish a paper you can cite if you found SBMLImporter helpful in your work.
