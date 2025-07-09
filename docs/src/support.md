# [Supported SBML Features](@id support)

SBMLImporter.jl supports a wide range of SBML features for models of level 2 or higher. This page provides an overview of the currently supported features, as well as the coverage of the SBML semantic and stochastic [test suites](https://github.com/sbmlteam/sbml-test-suite).

## Overview

In general, the following SBML features are supported:

* Events
* Rate rules
* Assignment rules
* Algebraic rules
* Dynamic compartment sizes
* Species and model conversion factors

Species can be defined in terms of concentration or amount, with units determined as follows:

* If `initialConcentration` is specified, the species is treated as a concentration.
* If `initialAmount` is specified, it is treated as an amount.
* If neither is set, and the species’ `substanceUnits` is `substance`, it defaults to amount.

The following SBML features are not currently supported:

* Delays (i.e., delay differential equations)
* Events with delay
* Events with priority
* Hierarchical (composite) models
* Fast reactions
* Parameter or species names that conflict with Julia constants (`pi`, `Inf`, `NaN`, `true`, `false`)
* Certain uncommon MathML expressions (e.g., `lt` with three arguments, `implies`, etc.)

## Semantic Test Suite Coverage

The [SBML semantic test suite](https://github.com/sbmlteam/sbml-test-suite/tree/release/cases/semantic) evaluates the correctness of deterministic ODE model simulations. SBMLImporter.jl currently passes 1266 out of 1821 test cases. A complete list of passed tests, and their corresponding `testTags` and `componentTags` can be found in [this](https://github.com/sebapersson/SBMLImporter.jl) CSV-file, which is automatically updated every time a new release of the package is made.

Overall, at least one test case passes for each of the following SBML `testTags`:

`Amount`, `NonUnityStoichiometry`, `BoundaryCondition`, `NonUnityCompartment`, `ConstantSpecies`, `LocalParameters`, `InitialValueReassigned`, `NonConstantParameter`, `0D-Compartment`, `NonConstantCompartment`, `MultiCompartment`, `ReversibleReaction`, `HasOnlySubstanceUnits`, `AssignedConstantStoichiometry`, `DefaultValue`, `Concentration`, `EventT0Firing`, `UncommonMathML`, `AssignedVariableStoichiometry`, `ConversionFactors`, `SpeciesReferenceInMath`, `VolumeConcentrationRates`, `NoMathML`, `L3v2MathML`, `BoolNumericSwap`.

And for each of the following SBML `componentTags`:

`Compartment`, `Species`, `Reaction`, `Parameter`, `FunctionDefinition`, `EventNoDelay`, `InitialAssignment`, `AssignmentRule`, `RateRule`, `AlgebraicRule`, `StoichiometryMath`, `CSymbolTime`, `CSymbolAvogadro`, `CSymbolRateOf`.

## Stochastic Test Suite Coverage

The [SBML stochastic test suite](https://github.com/sbmlteam/sbml-test-suite/tree/release/cases/stochastic) primarily evaluates the correctness of stochastic model simulations using an exact algorithm. SBMLImporter.jl currently passes 37 out of 100 test cases. The low pass rate is because cases 40–100 involve the SBML Distribution extension, which is currently not supported. As with the semantic test suite, a complete list of passed test cases along with their associated `testTags` and `componentTags` is available in [this CSV file](https://github.com/sebapersson/SBMLImporter.jl).

Overall, at least one test case passes for each of the following `testTags`:

`Amount`, `NonUnityStoichiometry`, `BoundaryCondition`, `NonUnityCompartment`, `ConstantSpecies`, `LocalParameters`, `InitialValueReassigned`.

And for each of the following SBML `componentTags`:

`Compartment`, `Species`, `Reaction`, `Parameter`, `FunctionDefinition`, `EventNoDelay`, `InitialAssignment`, `AssignmentRule`, `RateRule`, `AlgebraicRule`.

## Support for additional features

If SBMLImporter lacks support for a feature you would like to have, please file an issue on [GitHub](https://github.com/sebapersson/SBMLImporter.jl).
