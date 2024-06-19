function parse_rules!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for rule in libsbml_model.rules
        _parse_rule!(model_SBML, rule, libsbml_model)
    end
    return nothing
end

function _parse_rule!(model_SBML::ModelSBML, rule::SBML.AssignmentRule, libsbml_model::SBML.Model)::Nothing
    id, formula, have_ridents = _parse_rule_formula(rule, model_SBML, libsbml_model; assignment_rule=true)

    if _is_model_variable(id, model_SBML)
        variable = _get_model_variable(id, model_SBML)
        _add_rule_info!(variable, formula, have_ridents; assignment_rule=true)
        return nothing
    end

    # For StoichometryMath, assignment rules with variables names generatedId... are
    # created. These variables are only used in setting stoichometry for a reaction,
    # and do not correspond to specie, parameter and compartment.
    if occursin("generatedId", id)
        model_SBML.generated_ids[id] = formula
        # TODO : Can probably remove from list of assignment rules already here
        return nothing
    end

    # Per the standard, an assignment rule can create a new specie. Not recomended.
    @warn "Assignment rule creates new specie $(id). Happens when $(id) does " *
          "not correspond to any model specie, parameter, or compartment."
    c = first(keys(libsbml_model.compartments))
    model_SBML.species[id] = SpecieSBML(id, false, false, formula, formula, c, "", :Amount, false, true, false, false, have_ridents)
    return nothing
end
function _parse_rule!(model_SBML::ModelSBML, rule::SBML.RateRule, libsbml_model::SBML.Model)::Nothing
    id, formula, have_ridents = _parse_rule_formula(rule, model_SBML, libsbml_model; rate_rule=true)

    if _is_model_variable(id, model_SBML)
        variable = _get_model_variable(id, model_SBML)
        _add_rule_info!(variable, formula, have_ridents; rate_rule=true)
        return nothing
    end

    # Per the standard, a rate rule can create a new specie. Not recomended. If the rate
    # rule corresponds to a SpecieReference, the reference sets the initial value
    # TODO: Fix warning for SpecieReference (next round)
    @warn "Rate rate rule creates new specie $(id). Happens when $(id) does " *
          "not correspond to any model specie, parameter, or compartment."
    c = first(keys(libsbml_model.compartments))
    initial_value = _get_specieref_initial_value(id, libsbml_model)
    model_SBML.species[id] = SpecieSBML(id, false, false, initial_value, formula, c, "", :Amount, false, false, true, false, have_ridents)
    return nothing
end
function _parse_rule!(model_SBML::ModelSBML, rule::SBML.AlgebraicRule, libsbml_model::SBML.Model)::Nothing
    _, formula, _ = _parse_rule_formula(rule, model_SBML, libsbml_model; algebraic_rule=true)
    # As an algebraic rule gives the equation for a specie/parameter/compartment variable
    # is just a place-holder name
    if isempty(model_SBML.algebraic_rules)
        variable = "1"
    else
        variable = maximum(keys(model_SBML.algebraic_rules)) * "1"
    end
    model_SBML.algebraic_rules[variable] = "0 ~ " * formula
    return nothing
end

function _parse_rule_formula(rule::SBMLRule, model_SBML::ModelSBML, libsbml_model::SBML.Model; assignment_rule::Bool=false, rate_rule::Bool=false, algebraic_rule::Bool=false)::Tuple{String, String, Bool}
    @assert any([assignment_rule, algebraic_rule, rate_rule]) "No rule set to be parsed"

    variable = algebraic_rule ? "" : rule.variable
    rate_rule == true && push!(model_SBML.rate_rule_variables, variable)
    assignment_rule == true && push!(model_SBML.assignment_rule_variables, variable)

    math_expression = parse_math(rule.math, libsbml_model)
    have_ridents = _has_reactionid_ident(math_expression.math_idents, libsbml_model)

    formula = process_SBML_str_formula(math_expression.formula, model_SBML, libsbml_model; check_scaling=true, assignment_rule=assignment_rule, rate_rule=rate_rule, algebraic_rule=algebraic_rule, variable=variable)
    formula = isempty(formula) ? "0.0" : formula
    return variable, formula, have_ridents
end

function _add_rule_info!(specie::SpecieSBML, formula::String, have_ridents::Bool; assignment_rule::Bool=false, rate_rule::Bool=false)::Nothing
    specie.assignment_rule = assignment_rule
    specie.rate_rule = rate_rule
    if isempty(specie.initial_value) && !rate_rule
        specie.initial_value = formula
    end
    # If unit amount adjust for SBML equations per standard being given in conc.
    formula = _adjust_for_unit(formula, specie)
    specie.formula = formula
    specie.has_reaction_ids = have_ridents
    return nothing
end
function _add_rule_info!(variable::Union{ParameterSBML, CompartmentSBML}, formula::String, have_ridents::Bool; assignment_rule::Bool=false, rate_rule::Bool=false)::Nothing
    variable.assignment_rule = assignment_rule
    variable.rate_rule = rate_rule
    # For a rate-rule, the initial value is given by the old formula/value
    if rate_rule == true
        variable.initial_value = variable.formula
    else
        variable.initial_value = formula
    end
    variable.formula = formula
    variable.has_reaction_ids = have_ridents
    return nothing
end

function _get_specieref_initial_value(id::String, libsbml_model::SBML.Model)::String
    specieref = _which_specieref(id, libsbml_model)
    if isnothing(specieref)
        return "1.0"
    else
        return _parse_variable(specieref.stoichiometry; default="1")
    end
end

# TODO: This is a different beast, and comes later
function identify_algebraic_rule_variables!(model_SBML::ModelSBML)::Nothing
    # In case the model has algebraic rules some of the formulas (up to this point) are zero. To figure out
    # which variable check which might be applicable
    if isempty(model_SBML.algebraic_rules)
        return nothing
    end

    candidates = String[]
    for (specie_id, specie) in model_SBML.species
        if specie.rate_rule == true || specie.assignment_rule == true ||
           specie.constant == true
            continue
        end
        if specie_id ∈ model_SBML.species_in_reactions && specie.boundary_condition == false
            continue
        end
        if !(specie.formula == "0.0" || isempty(specie.formula))
            continue
        end

        # Check if specie-id actually occurs in any algebraic rule
        should_continue::Bool = false
        for (rule_id, rule) in model_SBML.algebraic_rules
            if replace_variable(rule, specie_id, "") != rule
                should_continue = false
            end
        end
        should_continue == true && continue

        push!(candidates, specie_id)
    end
    # To be set as algebraic rule variable a specie must not appear in reactions,
    # or be a boundary condition. Sometimes several species can fulfill this for
    # a single rule, in this case we choose the first valid
    if !isempty(candidates)
        model_SBML.species[candidates[1]].algebraic_rule = true
        push!(model_SBML.algebraic_rule_variables, candidates[1])
    end

    for (parameter_id, parameter) in model_SBML.parameters
        if parameter.rate_rule == true || parameter.assignment_rule == true ||
           parameter.constant == true
            continue
        end

        # Check if specie-id actually occurs in any algebraic rule
        should_continue::Bool = false
        for (rule_id, rule) in model_SBML.algebraic_rules
            if replace_variable(rule, parameter_id, "") != rule
                should_continue = false
            end
        end
        should_continue == true && continue

        parameter.algebraic_rule = true
        push!(model_SBML.algebraic_rule_variables, parameter.name)
    end

    for (compartment_id, compartment) in model_SBML.compartments
        if compartment.rate_rule == true || compartment.assignment_rule == true ||
           compartment.constant == true
            continue
        end

        # Check if specie-id actually occurs in any algebraic rule
        should_continue::Bool = false
        for (rule_id, rule) in model_SBML.algebraic_rules
            if replace_variable(rule, compartment_id, "") != rule
                should_continue = false
            end
        end
        should_continue == true && continue

        compartment.algebraic_rule = true
        push!(model_SBML.algebraic_rule_variables, compartment.name)
    end

    return nothing
end
