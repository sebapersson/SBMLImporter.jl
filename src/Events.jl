function parse_events!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    event_index = 1
    for (event_id, event) in libsbml_model.events
        trigger = _parse_trigger(event.trigger, model_SBML, libsbml_model)
        if isempty(trigger)
            continue
        end

        assignments = _parse_assignments(event.event_assignments, model_SBML, libsbml_model)

        if !isnothing(event_id)
            name = event_id
        else
            name = "event" * string(event_index)
            event_index += 1
        end
        model_SBML.events[name] = EventSBML(name, trigger, assignments, event.trigger.initial_value)
    end
    return nothing
end

# Rewrites triggers in events to the correct Julia syntax
function _parse_trigger(trigger_sbml::Union{Nothing, SBML.Trigger}, model_SBML::ModelSBML, libsbml_model::SBML.Model)::String
    math_expression = parse_math(trigger_sbml.math, libsbml_model)
    trigger = process_SBML_str_formula(math_expression.formula, model_SBML, libsbml_model; check_scaling=true, scale_rateof=false)
    trigger = replace_reactionid_formula(trigger, libsbml_model)

    # TODO : This info should live in math_expression for early exit
    # If trigger is a number event is triggered when value is not equal to zero
    is_number(trigger) && return trigger * " != 0"
    isempty(trigger) && return ""
    if trigger[1:2] == "if" || trigger[1:3] in ["and", "xor"]
        throw(SBMLSupport("Events with gated triggers (if, and, xor) are not supported"))
    end

    # If the compartment is an assignment rule it will be simplified away when calling
    # structurally_simplified, and thus cannot be found by the callback in ODE models.
    # Therefore, for these occurances rul formula is inserted
    trigger = _adjust_assignment_rule_variables(trigger, model_SBML)

    # TODO: This should not be needed, must fix when arrive at Callbacks.jl
    if occursin(r"<|>", trigger)
        trigger = replace(trigger, "<" => "≤")
        trigger = replace(trigger, ">" => "≥")
    end
    return trigger
end

function _parse_assignments(event_assignments::Vector{SBML.EventAssignment}, model_SBML::ModelSBML, libsbml_model::SBML.Model)::Vector{String}
    formulas, assigned_to = String[], String[]
    for event_assignment in event_assignments
        math_expression = parse_math(event_assignment.math, libsbml_model)
        isempty(math_expression.formula) && continue

        # For events t is accessed via the integrator interface
        formula = process_SBML_str_formula(math_expression.formula, model_SBML, libsbml_model; check_scaling=true)
        formula = replace_reactionid_formula(formula, libsbml_model)
        formula = replace_variable(formula, "t", "integrator.t")

        assign_to = event_assignment.variable
        if haskey(model_SBML.species, assign_to)
            # Formulas are given in concentration, so must adjust if assign_to is in amount
            formula = _adjust_for_unit(formula, model_SBML.species[assign_to])
        end

        if _is_model_variable(assign_to, model_SBML) == false
            @warn "Event creates new specie $(assign_to). Happens when $(assign_to) does " *
                  "not correspond to any model specie, parameter, or compartment."
            c = first(keys(libsbml_model.compartments))
            model_SBML.species[assign_to] = SpecieSBML(assign_to, false, false, "1.0", "", c, "", :Amount, false, false, false, false)
        end

        push!(assigned_to, assign_to)
        push!(formulas, formula)
    end

    # If event acts on a compartment, all concentrations in said compartment must
    # be adjusted via: conc_new = conc_old * V_old / V_new
    _adjust_event_compartment!!(formulas, assigned_to, model_SBML)

    # See comment in _parse_trigger for why needed
    for i in eachindex(formulas)
        formulas[i] = _adjust_assignment_rule_variables(formulas[i], model_SBML)
    end

    return assigned_to .* " = " .* formulas
end

function _adjust_event_compartment!!(formulas::Vector{String}, assigned_to::Vector{String}, model_SBML::ModelSBML)::Nothing
    if any([haskey(model_SBML.compartments, a) for a in assigned_to]) == false
        return nothing
    end

    i_compartments = findall(x -> haskey(model_SBML.compartments, x), assigned_to)
    for i in i_compartments
        compartment_id, compartment_formula = assigned_to[i], formulas[i]
        for (specie_id, specie) in model_SBML.species
            if specie.compartment != compartment_id || specie.unit == :Amount
                continue
            end
            if specie_id in assigned_to
                is = findfirst(x -> x == specie_id, assigned_to)
                formula = '(' * formulas[is] * ')'
            else
                is = nothing
                formula = specie_id
            end
            formula = formula * " * " * compartment_id * "/(" * compartment_formula * ')'
            if !isnothing(is)
                formulas[is] = formula
            else
                push!(assigned_to, specie_id)
                push!(formulas, formula)
            end
        end
    end
    return nothing
end
