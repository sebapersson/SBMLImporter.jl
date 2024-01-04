# Function generating callbacksets for time-depedent SBML piecewise expressions, as callbacks are more efficient than
# using ifelse (for example better integration stability, faster runtimes etc...)
function create_callbacks_SBML(system,
                               model_SBML::ModelSBML,
                               model_name::String)

    p_ode_problem_names::Vector{String} = string.(parameters(system))
    model_specie_names::Vector{String} = replace.(string.(states(system)), "(t)" => "")

    n_callbacks = length(keys(model_SBML.ifelse_parameters)) + length(keys(model_SBML.events))

    # Set function names
    model_name = replace(model_name, "-" => "_")
    callbacks = Vector{Any}(undef, n_callbacks)
    callback_str = ""

    # tstops function is added in the initialisation part of the first model callback (whether or not it is
    # a discrete or cont. callback)
    write_tstops, _tstops = create_tstops_function(model_SBML, model_specie_names, p_ode_problem_names::Vector{String})
    callback_str *= write_tstops

    k = 1
    # For ifelse parameter
    for parameter in keys(model_SBML.ifelse_parameters)
        first_callback = k == 1
        _affect, _cond, _callback, _active_t0 =  create_callback_ifelse(parameter, model_SBML,
                                                                        p_ode_problem_names, model_specie_names,
                                                                        first_callback, _tstops)
        callback_str *= _affect * _cond * _callback * _active_t0
        _affect_f = @RuntimeGeneratedFunction(Meta.parse(_affect))
        _cond_f = @RuntimeGeneratedFunction(Meta.parse(_cond))
        _active_t0_f = @RuntimeGeneratedFunction(Meta.parse(_active_t0))
        _get_cb = @RuntimeGeneratedFunction(Meta.parse(_callback))
        callbacks[k] = _get_cb(_cond_f, _affect_f, _active_t0_f)
        k += 1
    end

    # For classical SBML events
    for key in keys(model_SBML.events)
        first_callback = k == 1
        _affect, _cond, _callback, _initial_function = create_callback_SBML_event(key, model_SBML, p_ode_problem_names,
                                                                                  model_specie_names, first_callback, 
                                                                                  _tstops)
        callback_str *= _affect * _cond * _callback *  _initial_function

        _affect_f = @RuntimeGeneratedFunction(Meta.parse(_affect))

        # Condition can only be activated when going from false to true,
        # a variable checking this happens must be in the DiscreteCallback
        # condition
        if occursin("from_neg", _cond)
            __cond_f = @RuntimeGeneratedFunction(Meta.parse(_cond))
            _cond_f = let from_neg = [!model_SBML.events[key].trigger_initial_value]
                (u, t, integrator) -> __cond_f(u, t, integrator, from_neg)
            end
        else
            _cond_f = @RuntimeGeneratedFunction(Meta.parse(_cond))
        end

        # Some events can firse at time zero, if there is an _initial_function
        # ensure this
        _get_cb = @RuntimeGeneratedFunction(Meta.parse(_callback))
        if isempty(_initial_function)
            callbacks[k] = _get_cb(_cond_f, _affect_f)
        else
            _init_f = @RuntimeGeneratedFunction(Meta.parse(_initial_function))
            callbacks[k] = _get_cb(_cond_f, _affect_f, _init_f)
        end
        k += 1
    end

    if n_callbacks > 0
        _get_cbset = "function get_cbset(cbs)\n\treturn CallbackSet(" * prod("cbs[$i], " for i in 1:n_callbacks)[1:end-2] * ")\nend"
        get_cbset = @RuntimeGeneratedFunction(Meta.parse(_get_cbset))
        cbset = get_cbset(callbacks)
    else
        cbset = CallbackSet()
    end

    # Write callback to file if required, otherwise just return the string for the callback and tstops functions
    return cbset, callback_str
end


function create_callback_ifelse(parameter_name::String,
                                model_SBML::ModelSBML,
                                p_ode_problem_names::Vector{String},
                                model_specie_names::Vector{String},
                                first_callback::Bool, 
                                tstops::String)::Tuple{String, String, String, String}

    # Check if the event trigger depend on parameters which are to be i) estimated, or ii) if it depend on models state.
    # For i) we need to convert tspan. For ii) we cannot compute tstops (the event times) prior to starting to solve
    # the ODE so it most be cont. callback
    _condition, side_activated_with_time = model_SBML.ifelse_parameters[parameter_name]
    discrete_event = !(check_condition_has_states(_condition, model_specie_names))

    # Replace any state or parameter with their corresponding index in the ODE system to be comaptible with event
    # syntax
    for (i, specie_name) in pairs(model_specie_names)
        _condition = replace_variable(_condition, specie_name, "u["*string(i)*"]")
    end
    for (i, p_name) in pairs(p_ode_problem_names)
        _condition = replace_variable(_condition, p_name, "integrator.p["*string(i)*"]")
    end

    # Replace inequality with - (root finding cont. event) or with == in case of
    # discrete event
    replace_with = discrete_event == true ? "==" : "-"
    _condition_for_t0 = deepcopy(_condition) # Needed for checking active at t0 function
    _condition = replace(_condition, r"<=|>=|>|<" => replace_with)

    # Build the condition function
    condition_function = "\nfunction condition_" * parameter_name * "(u, t, integrator)\n"
    condition_function *= "\t" * _condition * "\nend\n"

    # Build the affect function
    i_ifelse_parameter = findfirst(x -> x == parameter_name, p_ode_problem_names)
    affect_function = "function affect_" * parameter_name * "!(integrator)\n"
    affect_function *= "\tintegrator.p[" * string(i_ifelse_parameter) * "] = 1.0\nend\n"

    # Build the callback formula
    callback_formula = "function get_callback_" * parameter_name * "(cond, affect!, init)"
    if discrete_event == false
        callback_formula *= "\tcb = ContinuousCallback(cond, affect!, "
    else
        callback_formula *= "\tcb = DiscreteCallback(cond, affect!, "
    end
    callback_formula *= "save_positions=(false, false), initialize=init)\n" # So we do not get problems with saveat in the ODE solver
    callback_formula *= "\treturn cb\nend"

    # Building a function which check if a callback is activated at time zero (as this is not something Julia will
    # check for us so must be done here). This is then passed
    side_inequality = side_activated_with_time == "right" ? "!" : "" # Check if true or false evaluates expression to true
    active_t0_function = "function is_active_t0_" * parameter_name * "!(c, u, t, integrator)\n"
    active_t0_function *= "\tt = 0.0 # Used to check conditions activated at t0=0\n" * "\tp[" * string(i_ifelse_parameter) * "] = 0.0 # Default to being off\n"
    active_t0_function *= "\tif " * side_inequality *"(" * condition_active_t0 * ")\n" * "\t\tintegrator.p[" * string(i_ifelse_parameter) * "] = 1.0\n\tend\n"
    if first_callback == false
        active_t0_function *= "end\n"
    else
        active_t0_function *= "\t_tstops = " * tstops * "\n"
        # Ensure tstops are within the integrator time-intervall
        active_t0_function *= "tstops=_tstops[@.((integrator.tdir * _tstops > integrator.tdir * integrator.sol.prob.tspan[1])*(integrator.tdir *_tstops < integrator.tdir * integrator.sol.prob.tspan[2]))]\n"
        active_t0_function *= "\ttstops = isempty(tstops) ? tstops : vcat(minimum(tstops) / 2.0, tstops)\n"
        active_t0_function *= "\tadd_tstop!.((integrator,), tstops)\nend\n"
    end

    return affect_function, condition_function, callback_formula, active_t0_function
end


function create_callback_SBML_event(event_name::String,
                                    model_SBML::ModelSBML,
                                    p_ode_problem_names::Vector{String},
                                    model_specie_names::Vector{String},
                                    first_callback::Bool, 
                                    tstops::String)::Tuple{String, String, String, String}

    event = model_SBML.events[event_name]
    _condition = event.trigger
    affects = event.formulas
    initial_value_cond = event.trigger_initial_value

    discrete_event = !(check_condition_has_states(_condition, model_specie_names))

    # If the event trigger does not contain a model state but fixed parameters it can at a maximum be triggered once.
    if discrete_event == false
        # If we have a trigger on the form a ≤ b then event should only be
        # activated when crossing the condition from left -> right. Reverse
        # holds for ≥
        affect_neg = occursin("≤", _condition)
    else
        # Build the SBML activation, which has a check to see that the condition crosses from false to
        # true, per SBML standard
        _condition = "\tcond = " * _condition * " && from_neg[1] == true\n\t\tfrom_neg[1] = !(" * _condition * ")\n\t\treturn cond"
    end

    # Replace any state or parameter with their corresponding index in the ODE system to be comaptible with event
    # syntax
    _condition_at_t0 = event.trigger
    for (i, specie_name) in pairs(model_specie_names)
        _condition = replace_variable(_condition, specie_name, "u["*string(i)*"]")
        _condition_at_t0 = replace_variable(_condition_at_t0, specie_name, "u["*string(i)*"]")
    end
    for (i, p_name) in pairs(p_ode_problem_names)
        _condition = replace_variable(_condition, p_name, "integrator.p["*string(i)*"]")
        _condition_at_t0 = replace_variable(_condition_at_t0, p_name, "integrator.p["*string(i)*"]")
    end
    # Build the condition function used in Julia file, for discrete checking that event indeed is coming from negative
    # direction
    if discrete_event == false
        _condition = replace(_condition, r"≤|≥" => "-")
        condition_function = "\nfunction condition_" * event_name * "(u, t, integrator)\n\t" * _condition * "\nend\n"

    elseif discrete_event == true
        condition_function = "\nfunction _condition_" * event_name * "(u, t, integrator, from_neg)\n"
        condition_function *= _condition * "\nend\n"
    end

    # Building the affect function (which can act on multiple states and/or parameters)
    affect_function = "function affect_" * event_name * "!(integrator)\n\tu_tmp = similar(integrator.u)\n\tu_tmp .= integrator.u\n"
    affect_function_body = "\tu_tmp = similar(integrator.u)\n\tu_tmp .= integrator.u\n"
    for (i, affect) in pairs(affects)
        # In RHS we use u_tmp to not let order affects, while in assigning LHS we use u
        affect_function1, affect_function2 = split(affect, "=")
        for j in eachindex(model_specie_names)
            affect_function1 = replace_variable(affect_function1, model_specie_names[j], "integrator.u["*string(j)*"]")
            affect_function2 = replace_variable(affect_function2, model_specie_names[j], "u_tmp["*string(j)*"]")
        end
        affect_function *= "\t" * affect_function1 * " = " * affect_function2 * '\n'
        affect_function_body *= "\t" * affect_function1 * " = " * affect_function2 * '\n' # For t0 events
    end
    # For Jump simulations to get correct SSA simulations
    affect_function *= "\tif integrator.sol.prob isa DiscreteProblem\n"
    affect_function *= "\t\treset_aggregated_jumps!(integrator)\n\tend\n"
    affect_function *= "end"
    for i in eachindex(p_ode_problem_names)
        affect_function = replace_variable(affect_function, p_ode_problem_names[i], "integrator.p["*string(i)*"]")
        affect_function_body = replace_variable(affect_function_body, p_ode_problem_names[i], "integrator.p["*string(i)*"]")
    end

    # In case the event can be activated at time zero build an initialisation function
    if discrete_event == true && initial_value_cond == false
        initial_value_str = "function init_" * event_name * "(c,u,t,integrator)\n"
        initial_value_str *= "\tcond = " * _condition_at_t0 * "\n"
        initial_value_str *= "\tif cond == true\n"
        #initial_value_str *= "\t" * affect_function_body * "\n\tend\n"
        initial_value_str *= "\t" * "" * "\n\tend\n"
    elseif discrete_event == false && initial_value_cond == false
        initial_value_str = "function init_" * event_name * "(c,u,t,integrator)\n"
        initial_value_str *= "\tcond = " * _condition_at_t0 * "\n" # We need a Bool not minus (-) condition
        initial_value_str *= "\tif cond == true\n"
        initial_value_str *= "\t" * affect_function_body * "\n\tend\n"
    else
        initial_value_str = ""
    end

    if first_callback == true
        if isempty(initial_value_str)
            initial_value_str = "function init_" * event_name * "(c,u,t,integrator)\n"
        end
        initial_value_str *= "\t_tstops = " * tstops * "\n"
        initial_value_str *= "tstops=_tstops[@.((integrator.tdir * _tstops > integrator.tdir * integrator.sol.prob.tspan[1])*(integrator.tdir *_tstops < integrator.tdir * integrator.sol.prob.tspan[2]))]\n"
        initial_value_str*= "\ttstops = isempty(tstops) ? tstops : vcat(minimum(tstops) / 2.0, tstops)\n"
        initial_value_str *= "\tadd_tstop!.((integrator,), tstops)\nend\n"
    elseif !isempty(initial_value_str)
        initial_value_str *= "end"
    end

    # Build the callback, consider initialisation if needed and direction for ContinuousCallback
    if isempty(initial_value_str)
        callback_formula = "function get_callback_" * event_name * "(cond, affect!)\n"
    else
        callback_formula = "function get_callback_" * event_name * "(cond, affect!, init)\n"
    end
    if discrete_event == false
        if affect_neg == true
            callback_formula *= "\tcb = ContinuousCallback(cond, nothing, affect!, "
        else
            callback_formula *= "\tcb = ContinuousCallback(cond, affect!, nothing, "
        end
        if !isempty(initial_value_str) 
            callback_formula *= "initialize=init, "
        end
    elseif discrete_event == true
        if isempty(initial_value_str)
            callback_formula *= "\tcb = DiscreteCallback(cond, affect!, "
        else
            callback_formula *= "\tcb = DiscreteCallback(cond, affect!, initialize=init, "
        end
    end
    callback_formula *= "save_positions=(false, false))\n" # So we do not get problems with saveat in the ODE solver
    callback_formula *= "\treturn cb\nend\n"

    return affect_function, condition_function, callback_formula, initial_value_str
end


# Function computing t-stops (time for events) for piecewise expressions using the symbolics package
# to symboically solve for where the condition is zero.
function create_tstops_function(model_SBML::ModelSBML,
                                model_specie_names::Vector{String},
                                p_ode_problem_names::Vector{String})::Tuple{String, String}

    write_tstops = "\nfunction __compute_tstops(integrator)\n"
    if isempty(model_SBML.ifelse_parameters) && isempty(model_SBML.events)
        write_tstops *= "\t return Float64[]\nend\n"
        return write_tstops, "Float64[]"
    end

    conditions = string.(vcat([model_SBML.ifelse_parameters[key][1] for key in keys(model_SBML.ifelse_parameters)], [e.trigger for e in values(model_SBML.events)]))
    tstops = Vector{String}(undef, length(conditions))
    tstops_to_float = Vector{String}(undef, length(conditions))
    for (i, condition) in pairs(conditions)

        # In case the activation formula contains a state we cannot precompute the t-stop time as it depends on
        # the actual ODE solution.
        if check_condition_has_states(condition, model_specie_names)
            tstops[i] = ""
            tstops_to_float[i] = ""
            continue
        end

        # We need to make the parameters and states symbolic in order to solve the condition expression
        # using the Symbolics package.
        _variables = "@variables t, "
        _variables *= prod(string.(collect(p_ode_problem_names)) .* ", " )[1:end-2] * " "
        _variables *= prod(string.(collect(model_specie_names)) .* ", " )[1:end-2]
        variables_symbolic = eval(Meta.parse(_variables))

        # Note - below order counts (e.g having < first results in ~= incase what actually stands is <=)
        _condition = replace(condition, r"≤|≥|<=|>=|<|>|==" => "~")
        condition_symbolic = eval(Meta.parse(_condition))

        # Expression for the time at which the condition is triggered
        local expression_time
        try
            expression_time = string.(Symbolics.solve_for(condition_symbolic, variables_symbolic[1], simplify=true))
        catch
            throw(SBMLSupport("Not possible to solve for time event is activated"))
        end

        for (i, specie_name) in pairs(model_specie_names)
            expression_time = replace_variable(expression_time, specie_name, "integrator.u["*string(i)*"]")
        end
        for (i, p_name) in pairs(p_ode_problem_names)
            expression_time = replace_variable(expression_time, p_name, "integrator.p["*string(i)*"]")
        end

        # dual_to_float is needed as tstops for the integrator cannot be of type Dual
        tstops[i] = expression_time # Used when we convert timespan
        i += 1
    end

    _tstops = "Float64[" * prod([isempty(_t) ? "" : _t * ", " for _t  in tstops])[1:end-2] * "]"
    write_tstops *= "\treturn " * _tstops * "\nend\n\n"

    return write_tstops, _tstops
end


function check_condition_has_states(condition::AbstractString, model_specie_names::Vector{String})::Bool
    for i in eachindex(model_specie_names)
        _condition = replace_variable(condition, model_specie_names[i], "")
        if _condition != condition
            return true
        end
    end
    return false
end
