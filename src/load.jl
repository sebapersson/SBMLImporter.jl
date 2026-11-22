"""
    load_SBML(path; massaction=false, kwargs...)

Import an [SBML](https://sbml.org/) model into a `ReactionNetwork` that can be simulated as
a `JumpProblem` for Gillespie simulations, a `SDEProblem` for chemical Langevin simulations,
or an `ODEProblem` for deterministic simulations.

The keyword `massaction` should be set to `true` if the model reactions follow chemical
[mass action](https://en.wikipedia.org/wiki/Law_of_mass_action) kinetics. This is because
the most efficient `JumpProblem` (Gillespie) simulators require mass action jumps, for more
details see [here](https://docs.sciml.ai/JumpProcesses/stable/jump_types/). For how to
check if a model follows mass action kinetics, see the FAQ in the documentation.

!!! note
    The `massaction` keyword argument is only relevant for jump simulations. If you import
    the model as a `SDEProblem` or an `ODEProblem`, you can (and should) ignore this
    keyword argument.

## Keyword arguments
- `complete::Bool=true`: Whether or not to mark the returned Catalyst `ReactionSystem` as
    complete. Only a complete system can be used for simulations, while only an incomplete
    system can be composed with other systems. For more details, see the Catalyst.jl
    [documentation](https://docs.sciml.ai/Catalyst/stable/).
- `ifelse_to_callback::Bool=true`: Rewrite `ifelse` (SBML piecewise) expressions to
  [callbacks](https://github.com/SciML/DiffEqCallbacks.jl). This improves
  simulation runtime and stability. Strongly recomended to set to `true`.
- `inline_assignment_rules::Bool=false`: Inline SBML assignment rules into model equations.
    Recommended for importing large models. **However**, if set to `true`, it is
    not possible to access assignment rule variables via `sol[:var]`.
- `inline_kineticlaw_parameters::Symbol=:inline`: How to handle SBML kinetic law
    parameters, which are parameters defined only in the kinetic law field of an SBML
    reaction. Can be one of:
    - `:inline` (default): Inline (hardcode) the parameter values into the model equations.
    - `:promote`: Promote the parameters to model parameters (with the same status as
        parameters in the SBML file), keeping their SBML ids. Parameters with the same id
        are merged if they have the same value, otherwise an error is thrown.
    - `:promote_qualified`: Promote all kinetic law parameters to model parameters with
        the reaction id as prefix (`<reaction id>_<parameter id>`), via the libSBML
        `promoteLocalParameters` converter. As the names are unique this also works when
        the same kinetic law parameter id is repeated with different values.
    Providing a `Bool` is deprecated (`true` corresponds to `:inline` and `false` to
    `:promote`).
- `write_to_file=false`: Write the parsed SBML model to a Julia file in the same
    directory as the SBML file.
- `model_as_string::Bool=false`: Whether the model (`path`) is provided as a string.
## Returns
- `rn`: A Catalyst reaction network/system that can be converted to a `JumpProblem`,
  `SDEProblem`, or `ODEProblem`. Initial values and parameter values can be accessed via
  Catalyst maps:
  - `get_u0_map(rn)` returns a map from species to initial values.
  - `get_parameter_map(rn)` returns a map from parameters to values.
- `cbset`: A `CallbackSet` containing SBML events and callbacks generated from SBML
  piecewise expressions (when `ifelse_to_callback=true`).
"""
function load_SBML(
        path::AbstractString; massaction::Bool = false, complete::Bool = true,
        ifelse_to_callback::Bool = true, write_to_file::Bool = false,
        inline_assignment_rules::Bool = false,
        inline_kineticlaw_parameters::Union{Bool, Symbol} = :inline,
        model_as_string::Bool = false
    )::Tuple{Catalyst.ReactionSystem, CallbackSet}
    model_SBML = parse_SBML(
        path, massaction; ifelse_to_callback = ifelse_to_callback,
        model_as_string = model_as_string,
        inline_assignment_rules = inline_assignment_rules,
        inline_kineticlaw_parameters = _get_kineticlaw_mode(inline_kineticlaw_parameters)
    )
    model_SBML_sys = _to_system_syntax(model_SBML, inline_assignment_rules, massaction)
    rn, specie_map, parameter_map = _get_reaction_system(model_SBML_sys, model_SBML)
    rn = Catalyst.set_u0_map(rn, specie_map)
    rn = Catalyst.set_parameter_map(rn, parameter_map)

    # Build callback functions
    cbset = create_callbacks(rn, model_SBML, model_SBML.name)

    if write_to_file == true && model_as_string == true
        @warn "If the model is provided as a string we do not support writing it to \
            file. Hence write_to_file argument is ignored"
    elseif write_to_file == true
        dirsave = _get_dir_save(write_to_file, model_as_string, path)
        write_reactionsystem(model_SBML_sys, dirsave, model_SBML)
    end

    if complete == true
        rn = Catalyst.complete(rn)
    end
    return rn, cbset
end

function _get_kineticlaw_mode(x::Union{Bool, Symbol})::Symbol
    if x isa Bool
        replacement = x ? ":inline" : ":promote"
        Base.depwarn(
            "Providing inline_kineticlaw_parameters as a Bool is deprecated, use the \
            Symbol $replacement instead of $x. Allowed values are :inline, :promote and \
            :promote_qualified.", :load_SBML; force = true
        )
        return x ? :inline : :promote
    end
    x in (:inline, :promote, :promote_qualified) && return x
    throw(
        ArgumentError(
            "Invalid inline_kineticlaw_parameters = :$x. Allowed values are :inline, \
            :promote and :promote_qualified"
        )
    )
end
