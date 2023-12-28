# [Supported SBML Features](@id support)

SBMLImporter supports many SBML features for models of level 2 or higher. Currently, excluding FBA models it successfully passes 1257 out of 1785 test cases. The failed test cases cover features currently not supported. Key features supported include:

- Events
- Rate rules
- Assignment rules
- Algebraic rules
- Dynamic compartment size
- Species and model conversion factors

Species can be specified in either concentration or amount. The unit determination is:

- If `initialConcentration` is set for a species, its unit will be set to concentration.
- If `initialAmount` is set for a species, it will be treated as being in amount.
- If neither is set, and the `substanceUnits` of the species is "substance," it is treated as being in amounts.

Currently SBMLImporter does not support the following:

* Models with empty SBML reactions
* Delay (creating a delay-differential-equations)
* Events with delay
* Events with priority
* Hierarchical models
* Fast reactions
* Parameter or species names corresponding to Julia constants (`pi`, `Inf`, `NaN`, `true`, `false`)
* Certain uncommon math expressions, such as `lt` with three arguments, `implies` etc...

Import might also fail for complicated nested piecewise expressions inside SBML functions.

## Support for additional features

If SBMLImporter lacks support for a feature you would like to have, please file an issue on GitHub.
