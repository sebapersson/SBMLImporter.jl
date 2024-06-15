function parse_rules!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for rule in libsbml_model.rules
        _parse_rule!(model_SBML, rule, libsbml_model)
    end
    return nothing
end

function _parse_rule!(model_SBML::ModelSBML, rule::SBML.AssignmentRule, libsbml_model::SBML.Model)::Nothing
    variable, formula = _parse_rule_formula(rule, model_SBML, libsbml_model; assignment_rule=true)

    if haskey(model_SBML.species, variable)
        _add_rule_info!(model_SBML.species[variable], formula; assignment_rule=true)
        return nothing
    end
    if haskey(model_SBML.parameters, variable)
        _add_rule_info!(model_SBML.parameters[variable], formula; assignment_rule=true)
        return nothing
    end
    if haskey(model_SBML.compartments, variable)
        _add_rule_info!(model_SBML.compartments[variable], formula; assignment_rule=true)
        return nothing
    end

    # For StoichometryMath, assignment rules with variables names generatedId... are
    # created. These variables are only used in setting stoichometry for a reaction,
    # and do not correspond to specie, parameter and compartment.
    if occursin("generatedId", variable)
        model_SBML.generated_ids[variable] = formula
        # TODO : Can probably remove from list of assignment rules already here
        return nothing
    end

    # Per the standard, an assignment rule can create a new specie. Not recomended.
    @warn "Assignment rule creates new specie $(variable). Happens when $(variable) does " *
          "not correspond to any model specie, parameter, or compartment."
    c = first(keys(libsbml_model.compartments))
    model_SBML.species[variable] = SpecieSBML(variable, false, false, formula, formula, c, "", :Amount, false, true, false, false)
    return nothing
end
function _parse_rule!(model_SBML::ModelSBML, rule::SBML.RateRule, libsbml_model::SBML.Model)::Nothing
    variable, formula = _parse_rule_formula(rule, model_SBML, libsbml_model; rate_rule=true)

    if haskey(model_SBML.species, variable)
        _add_rule_info!(model_SBML.species[variable], formula; rate_rule=true)
        return nothing
    end
    if haskey(model_SBML.parameters, variable)
        _add_rule_info!(model_SBML.parameters[variable], formula; rate_rule=true)
        return nothing
    end
    if haskey(model_SBML.compartments, variable)
        _add_rule_info!(model_SBML.compartments[variable], formula; rate_rule=true)
        return nothing
    end

    # Per the standard, a rate rule can create a new specie. Not recomended.
    @warn "Assignment rate rule creates new specie $(variable). Happens when $(variable) does " *
          "not correspond to any model specie, parameter, or compartment."
    c = first(keys(libsbml_model.compartments))
    model_SBML.species[variable] = SpecieSBML(variable, false, false, "1.0", formula, c, "", :Amount, false, false, true, false)
    return nothing
end
function _parse_rule!(model_SBML::ModelSBML, rule::SBML.AlgebraicRule, libsbml_model::SBML.Model)::Nothing
    _, formula = _parse_rule_formula(rule, model_SBML, libsbml_model; algebraic_rule=true)
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

function _parse_rule_formula(rule::SBMLRule, model_SBML::ModelSBML, libsbml_model::SBML.Model; assignment_rule::Bool=false, rate_rule::Bool=false, algebraic_rule::Bool=false)::Tuple{String, String}
    @assert any([assignment_rule, algebraic_rule, rate_rule]) "No rule set to be parsed"

    variable = algebraic_rule ? "" : rule.variable
    rate_rule == true && push!(model_SBML.rate_rule_variables, variable)
    assignment_rule == true && push!(model_SBML.assignment_rule_variables, variable)

    math_expression = parse_math(rule.math, libsbml_model)
    formula = replace_reactionid_formula(math_expression.formula, libsbml_model)
    formula = process_SBML_str_formula(formula, model_SBML, libsbml_model; check_scaling=true, assignment_rule=assignment_rule, rate_rule=rate_rule, algebraic_rule=algebraic_rule, variable=variable)
    formula = isempty(formula) ? "0.0" : formula
    return variable, formula
end

function _add_rule_info!(specie::SpecieSBML, formula::String; assignment_rule::Bool=false, rate_rule::Bool=false)::Nothing
    specie.assignment_rule = assignment_rule
    specie.rate_rule = rate_rule
    if isempty(specie.initial_value) && !rate_rule
        specie.initial_value = formula
    end
    # If unit amount adjust for SBML equations per standard being given in conc.
    if specie.unit == :Amount && specie.only_substance_units == false
        formula = "(" * formula * ")*" * specie.compartment
    end
    specie.formula = formula
    return nothing
end
function _add_rule_info!(variable::Union{ParameterSBML, CompartmentSBML}, formula::String; assignment_rule::Bool=false, rate_rule::Bool=false)::Nothing
    variable.assignment_rule = assignment_rule
    variable.rate_rule = rate_rule
    # For a rate-rule, the initial value is given by the old formula/value
    if rate_rule == true
        variable.initial_value = variable.formula
    else
        variable.initial_value = formula
    end
    variable.formula = formula
    return nothing
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
