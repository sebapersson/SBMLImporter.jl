# [FAQs](@id FAQ)

Certainly! Here's the revised **Markdown version** of the FAQ **without subheadings**, keeping it clean and linear while still clear and professional:

## [Why do I encounter installation problems?](@id install_fail)

SBMLImporter.jl is regularly tested on, and should be installable on Linux, macOS and Windows. If you encounter installation issues on these systems, we recommend checking the following two common causes:

1. an incorrectly installed or outdated Julia version,  
2. outdated package dependencies.

First, ensure that a supported Julia version is used. SBMLImporter is tested with Julia LTS version **1.10** and the latest stable version. Using an earlier version may result in installation failures. To reliably install and manage Julia versions across operating systems, we strongly recommend using [juliaup](https://github.com/JuliaLang/juliaup). If you are constrained to using an older Julia version, for example on an HPC cluster, and encounter problems, please file an issue on [GitHub](https://github.com/sebapersson/SBMLImporter.jl/issues).

Second, installation failures may result from outdated versions of SBMLImporter dependencies. For example, if SBMLImporter is installed into the global Julia environment, older versions of other packages may prevent the latest version from being installed. This can cause installation failures or break tutorials and example code. To avoid this, it is recommended to install SBMLImporter in a new, isolated environment. For example, to install it in an environment named `sbmlimporter_project`, run the following in the Julia REPL:

```julia
using Pkg
Pkg.activate("sbmlimporter_project")
Pkg.add("SBMLImporter")
# Add any additional packages as needed
```

If you need to install SBMLImporter into an existing environment and encounter issues, updating all packages may resolve the problem:

```julia
Pkg.update()
```

This is because SBMLImporter depends on numerous packages from the actively developed Julia [SciML ecosystem](https://sciml.ai/). New releases of these dependencies sometimes introduce breaking changes that are not always caught by test suites (e.g., see [this issue](https://github.com/SciML/Catalyst.jl/issues/1075)). In other words, SBMLImporter is not compatible with all versions of packages like Catalyst.jl, which can cause issues if an incompatible version is already installed in the environment.

Lastly, if you have tried everything above (e.g., you are using a recent Julia version and have updated your packages) and still experience installation issues, it is likely a bug in SBMLImporter. In this case, please open an issue on [GitHub](https://github.com/sebapersson/SBMLImporter.jl).

## Why do I get the error `Any[...] are either missing from the variable map or missing from the system's states/parameters list`?

This error probably occurs because the model has an SBML [assignment rule](https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/java-api/?org/sbml/libsbml/AssignmentRule.html). Assignment rules describe the assignment of an equation to a variable:

```julia
var = assignment # var ~ assignment in ModelingToolkit syntax
```

To be able to simulate the model assignment rules must be inlined into the model equations. This can be done with the `structural_simplify` function from [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl). For example, for an ODE model do:

```julia
using SBMLImporter, ModelingToolkit
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
    The `massaction` keyword argument is only relevant for jump simulations. If you import the model as a `SDEProblem` or an `ODEProblem`, you can (and should) ignore this keyword argument.

## Why does one of the SBML model parameters not appear among the system parameters?

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

If a parameter is set to have `constant="false"`, the importer must treat the parameter as a variable since it can change over time (explicitly depend on time), as Julia `ReactionSystem` parameters are assumed to time invariant. If the parameter should indeed be constant (e.g. at most be changed by an event), change the parameter in the model file:

```SBML
<parameter id="c" value="1.0" constant="true"/>
```

## How do I access SBML reaction ID and reaction name?

SBML reactions have both an `ID` and a `name` that can differ. When importing an SBML model, SBMLImporter stores these as `metadata` in every `Catalyst.Reaction`. This `metadata` can be accessed with the Catalyst `getmetadata` function. For example, to retrieve both the ID and name for the first reaction in a model, do:

```julia
using SBMLImporter, Catalyst
prn, cb = load_SBML(path_SBML)
sbml_reactions = reactions(prn.rn)
getmetadata(sbml_reactions[1], :id)
getmetadata(sbml_reactions[1], :name)
```

## How do I get the compartment of an SBML specie?

Every SBML `specie` has a `compartment`. SBMLImporter stores the compartment as `Catalyst.Species` metadata, which can be accessed via the `getcompartment` function. For more details, see the [API](@ref API).

## Why does my simulation fail with `DomainError` while the model imports fine?

This typically happens due to two reasons. Firstly, the model might contain a function call where the argument must be positive, such as:

```julia
log(specie)
```

Even though the model might be written such that `specie` should never go below 0 (e.g., the model follows mass action kinetics), numerical errors can cause `specie` to go below zero. Therefore, instead of encoding risky statements into the model such as `log(specie)`, it might be better to encode something like `log(specie + 1e-8)`.

Secondly `DomainError` might arise due to how SBMLImporter parses SBML piecewise expressions. Piecewise expressions are first parsed into `ifelse` functions:

```julia
ifelse(condition, x, y)
```

Which when `condition == true` evaluates to `x`, otherwise to `y`. In SBMLImporter `ifelse` expressions are further rewritten to callbacks (events). Hence, in the model equations the `ifelse` is converted to:

```julia
(1 - condition_bool) * x + condition_bool * y
```

Here `condition_bool` is assigned via a callback (event) to be `1` when `condition` is `true`, and to `0` when `condition` if `false`. This has the same functionality as an `ifelse`, but is numerically more stable because the integrator (e.g., ODE solver) actually stops at points `condition` changes. Instead, with `ifelse`, the integrator does not know an event happens and thus must typically reduce its step size to handle the sudden change in models dynamics. This increases simulation time, and reduces simulation stability (e.g. ODE solvers fail more frequently). However, sometimes the formulas in the `ifelse` might depend on the `condition`. For example, let's say we have:

```julia
ifelse(t > 1, 0, log(t - 1))
```

With the `ifelse` formulation, `log(t - 1)` is never evaluated until `t > 1`. With the callback formulation:

```julia
(1 - condition_bool) * 0 + condition_bool * log(t-1)
```

`log(t-1)` is evaluated for time-points `t < 1` which causes `DomainError`. This can be solved by simply not rewriting `ifelse` to callback when importing the model:

```julia
prn, cb = load_SBML(path; ifelse_to_callback = false)
```

If neither of the above solutions work, please file an issue on [GitHub](https://github.com/sebapersson/SBMLImporter.jl).

## Why did I get the `SBMLSupport` error?

You likely got this error because your SBML model contains an SBML feature that SBMLImporter does not support yet. An extensive list of supported features can be found [here](@ref support). If SBMLImporter lacks support for a feature you would like to have, please file an issue on [GitHub](https://github.com/sebapersson/SBMLImporter.jl).

## Do you support SBML export?

We currently do not support exporting Catalyst `ReactionSystem` to SBML. Most things are available in Julia for supporting SBML export though, so if anyone is interested, we would be happy to provide guidance for creating a SBMLExporter package.
