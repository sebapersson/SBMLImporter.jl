# SBMLImporter.jl
*Julia Importer for SBML models*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/dev/)
[![Build Status](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/sebapersson/SBMLImporter.jl/graph/badge.svg?token=1J7RKUPZOM)](https://codecov.io/gh/sebapersson/SBMLImporter.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

SBMLImporter is a Julia package for importing dynamic [Systems Biology Markup Language (SBML)](https://sbml.org/) models into either a `JumpProblem` for Gillespie simulations, a `SDEProblem` for chemical Langevin simulations, or an `ODEProblem` for deterministic simulations. Some major highlights of SBMLImporter are:

* It imports models into a [Catalyst.jl](https://github.com/SciML/Catalyst.jl) `ReactionSystem`. This allows for easy conversion to a `JumpProblem`, a `SDEProblem`, or an `ODEProblem`.
* It supports a majority of SBML features, such as dynamic compartments, events, rules, piecewise (ifelse) expressions, and units.
* It is thoroughly tested against both the SBML test suite and a large collection of published models.
* It integrates with [PEtab.jl](https://github.com/sebapersson/PEtab.jl) for fitting SBML models to data.

Additional information on features and tutorials can be found in the [documentation](https://sebapersson.github.io/SBMLImporter.jl/stable/).

## Installation

SBMLImporter is a registered Julia package and can be installed with the Julia package manager using:

```julia
julia> import Pkg; Pkg.add("SBMLImporter")
```

SBMLImporter is compatible with Julia 1.10 and above. For additional installation details, see the [documentation](https://sebapersson.github.io/SBMLImporter.jl/stable/).

## Other SBML related Julia packages

There are currently three other SBML related packages in Julia:

1. [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl) wraps a subset of the [libSBML](https://sbml.org/software/libsbml/) functionality and is used by SBMLImporter and other SBML related packages for parsing SBML models.
2. [COBREXA.jl](https://github.com/COBREXA/COBREXA.jl) is designed for constraint-based metabolic modeling. Constraint-based models, which are often referred to as flux-balance analysis (FBA) models, are not supported by SBMLImporter.
3. [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl) is the most similar package to SBMLImporter as it imports an SBML model into a `ReactionSystem`. Still, we recommend to use SBMLImporter as it supports more SBML features, has more efficient event handling, and supports efficient `JumpProblem` (Gillespie) simulations. An extensive list of differences can be found in the [documentation](https://sebapersson.github.io/SBMLImporter.jl/stable/differences/).

## Citation

We will soon publish a paper you can cite if you found SBMLImporter helpful in your work.
