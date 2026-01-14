# SBMLImporter.jl
*Julia Importer for SBML models*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sebapersson.github.io/SBMLImporter.jl/dev/)
[![Build Status](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sebapersson/SBMLImporter.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/sebapersson/SBMLImporter.jl/graph/badge.svg?token=1J7RKUPZOM)](https://codecov.io/gh/sebapersson/SBMLImporter.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

[Getting Started](https://sebapersson.github.io/SBMLImporter.jl/stable/tutorial/) |
[Documentation](https://sebapersson.github.io/SBMLImporter.jl/stable/) |
[Contributing](https://sebapersson.github.io/SBMLImporter.jl/stable/contributing/)

SBMLImporter is a Julia package for importing dynamic
[Systems Biology Markup Language(SBML)](https://sbml.org/) models into either a
`JumpProblem` for Gillespie simulations, a `SDEProblem` for chemical Langevin simulations,
or an `ODEProblem` for deterministic simulations. Major features are:

- Imports an SBML models into a [Catalyst.jl](https://github.com/SciML/Catalyst.jl)
  `ReactionSystem`. This allows for easy conversion to a `JumpProblem`, a `SDEProblem`, or
  an `ODEProblem`.
- Support for a majority of SBML features, such as dynamic compartments, events, rules,
  piecewise (`ifelse`) expressions, and units.
- Thoroughly tested against both the
  [SBML test suite](https://github.com/sbmlteam/sbml-test-suite) and a large collection of
  published models.
- Integrates with [PEtab.jl](https://github.com/sebapersson/PEtab.jl) for parameter
  estimating SBML models.

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

If you use SBMLImporter in work that is published, please cite the paper below:

```bibtex
@article{PEtabBioinformatics2025,
  title={PEtab.jl: advancing the efficiency and utility of dynamic modelling},
  author={Persson, Sebastian and Fr{\"o}hlich, Fabian and Grein, Stephan and Loman, Torkel and Ognissanti, Damiano and Hasselgren, Viktor and Hasenauer, Jan and Cvijovic, Marija},
  journal={Bioinformatics},
  volume={41},
  number={9},
  pages={btaf497},
  year={2025},
  publisher={Oxford University Press}
}
```
