"""
    SBML_to_ODESystem(path_SBML::AbstractString;
                      ifelse_to_callback::Bool=true,
                      write_to_file::Bool=false, 
                      verbose::Bool=true, 
                      return_all::Bool=false, 
                      model_as_string::Bool=false)

Parse an SBML model into a ModelingToolkit `ODESystem` and potentially convert events/piecewise to callbacks.

By default, `structurally_simplified` is called on the `ODESystem` before it is returned.

For information on simulating the `ODESystem`, refer to the documentation.

For testing `path_SBML` can be the model as a string if `model_as_string=true`

!!! note
    The number of returned arguments depends on whether the SBML model has events and/or piecewise expressions (see below).

## Arguments
- `path_SBML`: File path to a valid SBML file (level 2 or higher).
- `ifelse_to_callback=true`: Whether to rewrite `ifelse` (piecewise) expressions to callbacks; recommended for performance.
- `write_to_file=false`: Whether to write the parsed SBML model to a Julia file in the same directory as the SBML file.
- `verbose=true`: Whether or not to display information on the number of return arguments.
- `return_all=true`: Whether or not to return all possible arguments (see below), regardless of whether the model has events.
- `model_as_string=false` : Whether or not the model (`path_SBML`) is provided as str, mainly for testing.

## Returns
- `ode_system`: A ModelingToolkit `ODESystem` that can be converted into an `ODEProblem` and solved.
- `specie_map`: A species map setting initial values; together with the `ODESystem`, it can be converted into an `ODEProblem`.
- `parameter_map` A parameter map setting parameter values; together with the `ODESystem`, it can be converted into an `ODEProblem`.
- `cbset` - **only for models with events/piecewise expressions**: Callbackset (events) for the model.
- `get_tstops`- **Only for models with events/piecewise expressions**: Function computing time stops for discrete callbacks in the `cbset`.
- `ifelse_t0` - **Only for models with time-dependent ifelse (piecewise) expressions**: Functions checking and adjusting for callback-rewritten piecewise expressions that are active at `t=t0`.

## Examples
```julia
# Import and simulate model without events
using SBMLImporter
sys, specie_map, parameter_map = SBML_to_ODESystem(path_SBML)

using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Solve ODE with Rodas5P solver
sol = solve(prob, Rodas5P())
```
```julia
# Import a model with events
using SBMLImporter
sys, specie_map, parameter_map, cb, get_tstops = SBML_to_ODESystem(path_SBML)

using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Compute event times
tstops = get_tstops(prob.u0, prob.p)
sol = solve(prob, Rodas5P(), tstops=tstops, callback=callbacks)
```
"""                 
function SBML_to_ODESystem(path_SBML::T;
                           ifelse_to_callback::Bool=true,
                           write_to_file::Bool=false, 
                           verbose::Bool=true, 
                           ret_all::Bool=false, 
                           model_as_string::Bool=false) where T <: AbstractString

    # Intermediate model representation of a SBML model which can be processed into 
    # an ODESystem                           
    model_SBML = build_SBML_model(path_SBML; ifelse_to_callback=ifelse_to_callback, model_as_string=model_as_string)

    # If model is written to file save it in the same directory as the SBML-file
    if model_as_string == false
        dir_save = joinpath(splitdir(path_SBML)[1], "SBML")
    else
        dir_save = joinpath(@__DIR__, "SBML")
    end
    if write_to_file == true && !isdir(dir_save)
        mkdir(dir_save)
    end
    path_save_model = joinpath(dir_save, model_SBML.name * ".jl")

    # Build the ODESystem
    _model = odesystem_from_SBML(model_SBML, path_save_model, model_SBML.name, write_to_file)
    _get_ode_system = @RuntimeGeneratedFunction(Meta.parse(_model))
    _ode_system, specie_map, parameter_map = _get_ode_system("https://xkcd.com/303/") # Argument needed by @RuntimeGeneratedFunction
    # For a DAE we need special processing to rewrite it into an ODE 
    if isempty(model_SBML.algebraic_rules)
        ode_system = structural_simplify(_ode_system)
    else
        ode_system = structural_simplify(dae_index_lowering(_ode_system))
    end

    # Build callback functions 
    _callbacks, _tstops = create_callbacks_SBML(ode_system, model_SBML, model_SBML.name)
    get_callbacks = @RuntimeGeneratedFunction(Meta.parse(_callbacks))
    cbset, ifelse_t0 = get_callbacks("https://xkcd.com/2694/") # Argument needed by @RuntimeGeneratedFunction
    compute_tstops = @RuntimeGeneratedFunction(Meta.parse(_tstops))

    # if model is written to file write the callback
    if write_to_file == true
        path_save = joinpath(dir_save, model_SBML.name * "_callbacks.jl")
        io = open(path_save, "w")
        write(io, _callbacks * "\n\n")
        write(io, _tstops)
        close(io)
    end

    if ret_all == true
        return ode_system, specie_map, parameter_map, cbset, compute_tstops, ifelse_t0
    end

    # Depending on model return what is needed to perform forward simulations 
    if isempty(model_SBML.events) && isempty(model_SBML.ifelse_bool_expressions)
        return ode_system, specie_map, parameter_map
    end 
    if isempty(model_SBML.ifelse_bool_expressions) && !isempty(model_SBML.events)
        verbose && @info "SBML model with events - output returned as odesys, specie_map, parameter_map, cbset, compute_tstops\nFor how to simulate model see documentation"
        return ode_system, specie_map, parameter_map, cbset, compute_tstops
    end
    if isempty(model_SBML.ifelse_parameters) && !isempty(cbset)
        verbose && @info "SBML model with piecewise rewritten to event - output returned as odesys, specie_map, parameter_map, cbset, compute_tstops, cb_t0\nFor how to simulate model see documentation"
        return ode_system, specie_map, parameter_map, cbset, compute_tstops, ifelse_t0
    end
end


function odesystem_from_SBML(model_SBML::ModelSBML,
                             path_save_model::String,
                             model_name::String,
                             write_to_file::Bool)::String

    _eqs_write = "\teqs = [\n"
    _ODESystem_write = "\t@named sys = ODESystem(eqs, t, species, parameters)"

    # Check if model is empty of derivatives if the case add dummy state to be able to
    # simulate the model
    if ((isempty(model_SBML.species) || sum([!s.assignment_rule for s in values(model_SBML.species)]) == 0) &&
        (isempty(model_SBML.parameters) || sum([p.rate_rule for p in values(model_SBML.parameters)]) == 0) &&
        (isempty(model_SBML.compartments) || sum([c.rate_rule for c in values(model_SBML.compartments)]) == 0))

        model_SBML.species["foo"] = SpecieSBML("foo", false, false, "1.0", "0.0", "1.0", "", :Amount,
                                                  false, false, false, false)
    end

    _species_write, _species_write_array, _specie_map_write = get_specie_map(model_SBML)
    _parameters_write, _parameters_write_array, _parameter_map_write = get_parameter_map(model_SBML)

    #=
        Build the model equations
    =#
    # Species
    for (specie_id, specie) in model_SBML.species

        if specie.algebraic_rule == true
            continue
        end

        formula = isempty(specie.formula) ? "0.0" : specie.formula
        if specie.assignment_rule == true
            eq = specie_id * " ~ " * formula
        else
            eq = "D(" * specie_id * ") ~ " * formula
        end
        _eqs_write *= "\t" * eq * ",\n"
    end
    # Parameters
    for (parameter_id, parameter) in model_SBML.parameters

        if parameter.constant == true || parameter.algebraic_rule == true
            continue
        end

        if parameter.rate_rule == false
            eq = parameter_id * " ~ " * parameter.formula
        else
            eq = "D(" * parameter_id * ") ~ " * parameter.formula
        end
        _eqs_write *= "\t" * eq * ",\n"
    end
    # Compartments
    for (compartment_id, compartment) in model_SBML.compartments

        if compartment.constant == true || compartment.algebraic_rule == true
            continue
        end

        if compartment.rate_rule == false
            eq = compartment_id * " ~ " * compartment.formula
        else
            eq = "D(" * compartment_id * ") ~ " * compartment.formula
        end
        _eqs_write *= "\t" * eq * ",\n"
    end
    # Algebraic rules
    for rule_formula in values(model_SBML.algebraic_rules)
        _eqs_write *= "\t" * rule_formula * ",\n"
    end
    _eqs_write *= "\t]"

    ### Writing to file
    model_name = replace(model_name, "-" => "_")
    io = IOBuffer()
    println(io, "function get_ODESystem_" * model_name * "(foo)")
    println(io, "\t# Model name: " * model_name)
    println(io, "")

    println(io, _species_write)
    println(io, _species_write_array)
    println(io, "")
    println(io, _parameters_write)
    println(io, _parameters_write_array)
    println(io, "")
    println(io, "    D = Differential(t)")
    println(io, "")
    println(io, _eqs_write)
    println(io, "")
    println(io, _ODESystem_write)
    println(io, "")
    println(io, _specie_map_write)
    println(io, "")
    println(io, "\t# SBML file parameter values")
    println(io, _parameter_map_write)
    println(io, "")
    println(io, "    return sys, specie_map, parameter_map")
    println(io, "")
    println(io, "end")
    model_str = String(take!(io))
    close(io)

    # In case user request file to be written
    if write_to_file == true
        open(path_save_model, "w") do f
            write(f, model_str)
        end
    end
    return model_str
end
