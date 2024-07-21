function create_callbacks(system, model_SBML::ModelSBML, model_name::String;
                          convert_stpan::Bool = false)::CallbackSet
    model_name = replace(model_name, "-" => "_")
    specie_ids = replace.(string.(unknowns(system)), "(t)" => "")
    # If parameters(system) is empty Any[] is returned with string.(...)
    parameter_ids = string.(parameters(system))
    parameter_ids = parameter_ids == Any[] ? String[] : parameter_ids

    n_callbacks = length(keys(model_SBML.events))
    callbacks = Vector{Union{DiscreteCallback, ContinuousCallback}}(undef, n_callbacks)

    # To reduce the number of returned functions by load_SBML tstops for any potential
    # DiscreteCallback are computed via a get_tstops functions, which is executed in the
    # initialisation part of the first callback
    tstops = _get_tstops(model_SBML, specie_ids, parameter_ids)

    first_callback::Bool, k = true, 1
    for event in values(model_SBML.events)
        discrete_callback = !(_condition_has_species(event.trigger, specie_ids))
        condition_f = _get_callback_condition(event, specie_ids, parameter_ids)
        affect_f! = _get_callback_affect(event, specie_ids, parameter_ids)
        init_f!, has_init_f = _get_callback_init(event, specie_ids, parameter_ids, tstops,
                                                 first_callback)
        event_direction = _get_event_direction(event, discrete_callback)

        if has_init_f == false
            _init_f! = SciMLBase.INITIALIZE_DEFAULT
        else
            _init_f! = init_f!
        end
        if discrete_callback == true
            callbacks[k] = DiscreteCallback(condition_f, affect_f!, initialize = _init_f!,
                                            save_positions = (false, false))
        elseif discrete_callback == false && event_direction == :left
            callbacks[k] = ContinuousCallback(condition_f, nothing, affect_f!,
                                              initialize = _init_f!,
                                              save_positions = (false, false))
        elseif discrete_callback == false && event_direction == :right
            callbacks[k] = ContinuousCallback(condition_f, affect_f!, nothing,
                                              initialize = _init_f!,
                                              save_positions = (false, false))
        else
            callbacks[k] = ContinuousCallback(condition_f, affect_f!, initialize = _init_f!,
                                              save_positions = (false, false))
        end
        k += 1
        first_callback = false
    end
    return CallbackSet(callbacks...)
end

function _get_callback_condition(event::EventSBML, specie_ids::Vector{String},
                                 parameter_ids::Vector{String})::Function
    condition = event.trigger
    discrete_callback = !(_condition_has_species(condition, specie_ids))

    # If the event is a rewritten ifelse the trigger operator must be rewritten to == or -
    # to ensure the event is properly triggered
    if event.is_ifelse == true
        if discrete_callback == true
            condition = replace(condition, r"<=|>=|>|<|≤|≥" => "==")
        else
            condition = replace(condition, r"<=|>=|>|<|≤|≥" => "-")
        end
    end
    condition = _ids_to_callback_syntax(condition, specie_ids, :specie; integrator = false)
    condition = _ids_to_callback_syntax(condition, parameter_ids, :parameter)

    condition_call = _template_condition(condition, discrete_callback, event.name)

    # Condition can only be activated when going from false to true, therefore an extra
    # captured Bool-variable for DiscreteCallback is added.
    if discrete_callback == true
        _condition_f = @RuntimeGeneratedFunction(Meta.parse(condition_call))
        condition_f = let from_neg = [!event.trigger_initial_value]
            (u, t, integrator) -> _condition_f(u, t, integrator, from_neg)
        end
    else
        condition_f = @RuntimeGeneratedFunction(Meta.parse(condition_call))
    end
    return condition_f
end

function _get_callback_affect(event::EventSBML, specie_ids::Vector{String},
                              parameter_ids::Vector{String})::Function
    affect_call = _template_affect(event, specie_ids, parameter_ids)
    return @RuntimeGeneratedFunction(Meta.parse(affect_call))
end

function _get_callback_init(event::EventSBML, specie_ids::Vector{String},
                            parameter_ids::Vector{String}, tstops::String,
                            first_callback::Bool)::Tuple{Function, Bool}
    discrete_callback = !(_condition_has_species(event.trigger, specie_ids))
    condition = event.trigger
    condition = _ids_to_callback_syntax(condition, specie_ids, :specie; integrator = false)
    condition = _ids_to_callback_syntax(condition, parameter_ids, :parameter)

    affect_body = _template_affect(event, specie_ids, parameter_ids; only_body = true)

    init_call = _template_init(event, condition, affect_body, tstops, first_callback,
                               discrete_callback)
    init_f! = @RuntimeGeneratedFunction(Meta.parse(init_call))
    has_init = event.trigger_initial_value == false || first_callback
    return init_f!, has_init
end

function _get_event_direction(event::EventSBML, discrete_callback::Bool)::Symbol
    # For continious callback if we have a trigger on the form a ≤ b then event should only
    # be activated when crossing the condition from left -> right (a gets bgger than b).
    # Reverse holds for ≥. If the trigger has == or - the event is activated from both
    # directions.
    if discrete_callback == false && event.is_ifelse == false
        if any(occursin.(["<=", "<", "≤"], event.trigger))
            event_direction = :left
        elseif any(occursin.([">=", ">", "≥"], event.trigger))
            event_direction = :right
        else
            event_direction = :none
        end
    else
        event_direction = :none
    end
    return event_direction
end

function _get_tstops(model_SBML::ModelSBML, specie_ids::Vector{String},
                     parameter_ids::Vector{String})::String
    tstops = String[]
    for event in values(model_SBML.events)
        if _condition_has_species(event.trigger, specie_ids)
            continue
        end

        # In this case we can solve for the event time via Symbolics
        variables = "ModelingToolkit.@variables t, "
        variables *= prod(parameter_ids .* ", ")[1:(end - 2)] * " "
        variables *= prod(specie_ids .* ", ")[1:(end - 2)] * " "
        variables_symbolic = eval(Meta.parse(variables))
        # In Symbolics syntax "~" gives equality
        condition = replace(event.trigger, r"≤|≥|<=|>=|<|>|==" => "~")
        condition_symbolic = eval(Meta.parse(condition))
        local tstop
        try
            tstop = string.(Symbolics.solve_for(condition_symbolic, variables_symbolic[1],
                                                simplify = true))
        catch
            throw(SBMLSupport("Not possible to solve for time event $(event.name) is activated"))
        end

        # In events species and parameters are given via integrator.u and integrator.p
        tstop = _ids_to_callback_syntax(tstop, specie_ids, :specie)
        tstop = _ids_to_callback_syntax(tstop, parameter_ids, :parameter)
        push!(tstops, tstop)
    end
    return _template_tstops(tstops)
end

function _condition_has_species(condition::String, specie_ids::Vector{String})::Bool
    for specie_id in specie_ids
        if condition != _replace_variable(condition, specie_id, "")
            return true
        end
    end
    return false
end

# TODO: When work also use interface for species in assignment
function _ids_to_callback_syntax(formula::String, ids::Vector{String}, id_type::Symbol;
                                 integrator::Bool = true, utmp::Bool = false)::String
    @assert id_type in [:specie, :parameter] "Invalid id type $id_type when parsing to callback syntax"
    if id_type == :parameter
        @assert integrator == true "Parameter not given via integrator in callback"
    end
    for (i, id) in pairs(ids)
        if id_type == :specie
            replace_with = utmp ? "utmp[" * string(i) * "]" : "u[" * string(i) * "]"
        elseif id_type == :parameter
            replace_with = "ps[:" * id * "]"
        end
        if integrator == true
            replace_with = "integrator." * replace_with
        end
        formula = _replace_variable(formula, id, replace_with)
    end
    return formula
end
