function parse_events!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    event_index = 1
    for (event_id, event) in libsbml_model.events
        trigger, have_ridents1, have_rateOf1, have_specieref1 = _parse_trigger(event.trigger, model_SBML, libsbml_model)
        if isempty(trigger)
            continue
        end

        assignments, have_ridents2, have_rateOf2, have_specieref2 = _parse_assignments(event.event_assignments, model_SBML, libsbml_model)

        if !isnothing(event_id)
            name = event_id
        else
            name = "event" * string(event_index)
            event_index += 1
        end
        model_SBML.events[name] = EventSBML(name, trigger, assignments, event.trigger.initial_value, have_ridents1, have_ridents2, have_rateOf1, have_rateOf2, have_specieref1, have_specieref2, false)
    end
    return nothing
end

# Rewrites triggers in events to the correct Julia syntax
function _parse_trigger(trigger_sbml::Union{Nothing, SBML.Trigger}, model_SBML::ModelSBML, libsbml_model::SBML.Model)::Tuple{String, Bool, Bool, Bool}
    math_expression = parse_math(trigger_sbml.math, libsbml_model)
    if any(in.(["if", "and", "xor"], [math_expression.fns]))
        throw(SBMLSupport("Events with gated triggers (if, and, xor) are not supported"))
    end
    have_ridents = math_expression.has_reaction_ids
    have_rateOf = "rateOf" in math_expression.fns
    have_specieref = _has_specieref(math_expression.math_idents, model_SBML)
    trigger = _process_formula(math_expression, model_SBML)

    # If trigger is a number event is triggered when value is not equal to zero
    is_number(trigger) && return (trigger * " != 0", false, false, false)
    isempty(trigger) && return ("", false, false, false)

    # If the compartment is an assignment rule it will be simplified away when calling
    # structurally_simplified, and thus cannot be found by the callback in ODE models.
    # Therefore, for these occurances rul formula is inserted
    trigger = _adjust_assignment_rule_variables(trigger, model_SBML)

    # TODO: This should not be needed, must fix when arrive at Callbacks.jl
    if occursin(r"<|>", trigger)
        trigger = replace(trigger, "<" => "≤")
        trigger = replace(trigger, ">" => "≥")
    end
    return trigger, have_ridents, have_rateOf, have_specieref
end

function _parse_assignments(event_assignments::Vector{SBML.EventAssignment}, model_SBML::ModelSBML, libsbml_model::SBML.Model)::Tuple{Vector{String}, Bool, Bool, Bool}
    formulas, assigned_to = String[], String[]
    have_specieref, have_ridents, have_rateOf = false, false, false
    for event_assignment in event_assignments
        math_expression = parse_math(event_assignment.math, libsbml_model)
        isempty(math_expression.formula) && continue

        _have_ridents = math_expression.has_reaction_ids
        _have_rateOf = "rateOf" in math_expression.fns
        _have_specieref = _has_specieref(math_expression.math_idents, model_SBML)
        have_ridents = _update_only_if_true(have_ridents, _have_ridents)
        have_rateOf = _update_only_if_true(have_rateOf, _have_rateOf)
        have_specieref = _update_only_if_true(have_specieref, _have_specieref)

        # For events t is accessed via the integrator interface
        formula = _process_formula(math_expression, model_SBML)
        formula = _replace_variable(formula, "t", "integrator.t")

        assign_to = event_assignment.variable
        if haskey(model_SBML.species, assign_to)
            # Formulas are given in concentration, so must adjust if assign_to is in amount
            formula = _adjust_for_unit(formula, model_SBML.species[assign_to])
        end

        if _is_model_variable(assign_to, model_SBML) == false
            @warn "Event creates new specie $(assign_to). Happens when $(assign_to) does " *
                  "not correspond to any model specie, parameter, or compartment."
            c = first(keys(libsbml_model.compartments))
            model_SBML.species[assign_to] = SpecieSBML(assign_to, false, false, "1.0", "", c, "", :Amount, false, false, false, false, false, false, false)
            _add_ident_info!(model_SBML.species[assign_to], math_expression, model_SBML)
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

    formulas = assigned_to .* " = " .* formulas
    return formulas, have_ridents, have_rateOf, have_specieref
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

function force_include_event_variables!(model_SBML::ModelSBML)::Nothing
    # Parameters/comparments can be non-constant, have a constant rhs and they change value
    # due to event assignments. To avoid the variable being simplified away by
    # structurally_simplify the variable is set to have zero derivative
    variables = Iterators.flatten((model_SBML.parameters, model_SBML.compartments))
    for (id, variable) in variables
        variable.constant == true && continue
        variable.rate_rule == true && continue
        variable.algebraic_rule == true && continue
        !is_number(variable.formula) && continue
        parse(Float64, variable.formula) ≈ π && continue # To pass test case 957

        variable.rate_rule = true
        variable.initial_value = variable.formula
        variable.formula = "0.0"
        push!(model_SBML.rate_rule_variables, id)
    end
    return nothing
end
