# [Supported SBML Features](@id support)

SBMLImporter provides extensive support for Systems Biology Markup Language (SBML) features. It successfully passes ADD test cases (excluding hierarchical models and FBA models). Key features supported include:

- Events
- Rate
- Assignment
- Algebraic rules
- Dynamic compartment size
- Species conversion factors

Species can be specified in either concentration or amount. The unit determination is:

- If `initialConcentration` is set for a species, its unit will be set to concentration.
- If `initialAmount` is set for a species, it will be treated as being in amount.
- If neither is set, and the `substanceUnits` of the species is "substance," it is treated as being in amounts.

Currently SBMLImporter does not support the following:

- Events with priority
- Events with delay
- Delay (creating a delay-differential-equation)
- Fast reactions
- Parameter or species names corresponding to Julia constants (`pi`, `NaN`, `true`, `false`)

Import might also fail for complicated nested piecewise expressions inside functions.

## Support for additional features

If SBMLImporter lacks support for a feature you would like, please file an issue on GitHub.