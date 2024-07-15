# [FAQs](@id FAQ)

## Why do I get the error `Any[...] are either missing from the variable map or missing from the system's states/parameters list`?

This error probably occurs because the model has an SBML [assignment rule](https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/java-api/?org/sbml/libsbml/AssignmentRule.html). Assignment rules describe the assignment of an equation to a variable:

```julia
var = assignment # var ~ assignment in ModelingToolkit syntax
```

To be able to simulate the model assignment rules must be inlined into the model equations. This can be done with the `structural_simplify` function from [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl). For example, for an ODE model do:

```julia
using ModelingToolkit
prn, cb = load_SBML(path)
sys = structural_simplify(convert(ODESystem, prn.rn))
oprob = ODEProblem(sys, prn.u0, tspan, prn.p)
```

## How can I check if my model follows mass action kinetics?

For efficient jump simulations (Gillespie type), the model should ideally follow chemical mass action kinetics. To inform the importer about mass action kinetics simply set the keyword argument `massaction=true` in `load_SBML`. Now, if you are unsure whether the model follows mass action kinetics, simply provide `massaction=true`. If the model does not adhere to mass action kinetics, a warning is thrown:

```
Warning: That the system is massaction was provided, however, the reaction r is not massaction. It is parsed as non massaction reaction, which can negatively impact simulation times
```

!!! note
    The `massaction` argument is only relevant for jump simulations. If you import the model as a `SDEProblem` or an `ODEProblem`, you can (and should) ignore this keyword argument.

## Why does one of the model parameters not appear among the system parameters?

If one of the model parameters does not appear among the system parameters:

```julia
prn, cb = load_SBML(path)
parameters(prn)
```

but appears in the model variables:

```julia
unknowns(prn)
```

this is likely because the parameter is not set as constant in the SBML file, e.g.:

```SBML
<parameter id="c" value="1.0" constant="false"/>
```

If a parameter is set to have `constant="false"`, the importer must treat the parameter as a variable since it can change over time (explicitly depend on time), because `ReactionSystem` parameters are assumed to time invariant. If the parameter should indeed be constant, change the parameter in the model file:

```SBML
<parameter id="c" value="1.0" constant="true"/>
```

## Why did I get the `SBMLSupport` error?

You likely got this error because your SBML model contains an SBML feature that SBMLImporter does not support yet. An extensive list of supported features can be found [here](@ref support). If SBMLImporter lacks support for a feature you would like to have, please file an issue on [GitHub](https://github.com/sebapersson/SBMLImporter.jl).

## Do you support SBML export?

We currently do not support exporting Catalyst `ReactionSystem` to SBML. Most things are available in Julia for supporting SBML export though, so if anyone is interested, we would be happy to provide guidance for creating a SBMLExporter package.
