function parse_rules!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for rule in libsbml_model.rules
        _parse_rule!(model_SBML, rule, libsbml_model)
    end
    return nothing
end

function _parse_rule!(model_SBML::ModelSBML, rule::SBML.AssignmentRule, libsbml_model::SBML.Model)::Nothing
    id, formula, math_expression = _parse_rule_formula(rule, model_SBML, libsbml_model; assignment_rule=true)
    have_ridents = _has_reactionid_ident(math_expression.math_idents, libsbml_model)
    have_rateOf = "rateOf" in math_expression.fns

    if _is_model_variable(id, model_SBML)
        variable = _get_model_variable(id, model_SBML)
        _add_rule_info!(variable, formula, have_ridents, have_rateOf; assignment_rule=true)
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
    model_SBML.species[id] = SpecieSBML(id, false, false, formula, formula, c, "", :Amount, false, true, false, false, have_ridents, have_rateOf)
    return nothing
end
function _parse_rule!(model_SBML::ModelSBML, rule::SBML.RateRule, libsbml_model::SBML.Model)::Nothing
    id, formula, math_expression = _parse_rule_formula(rule, model_SBML, libsbml_model; rate_rule=true)
    have_ridents = _has_reactionid_ident(math_expression.math_idents, libsbml_model)
    have_rateOf = "rateOf" in math_expression.fns

    if _is_model_variable(id, model_SBML)
        variable = _get_model_variable(id, model_SBML)
        _add_rule_info!(variable, formula, have_ridents, have_rateOf; rate_rule=true)
        return nothing
    end

    # Per the standard, a rate rule can create a new specie. Not recomended. If the rate
    # rule corresponds to a SpecieReference, the reference sets the initial value
    # TODO: Fix warning for SpecieReference (next round)
    @warn "Rate rate rule creates new specie $(id). Happens when $(id) does " *
          "not correspond to any model specie, parameter, or compartment."
    c = first(keys(libsbml_model.compartments))
    initial_value = _get_specieref_initial_value(id, libsbml_model)
    model_SBML.species[id] = SpecieSBML(id, false, false, initial_value, formula, c, "", :Amount, false, false, true, false, have_ridents, have_rateOf)
    return nothing
end
function _parse_rule!(model_SBML::ModelSBML, rule::SBML.AlgebraicRule, libsbml_model::SBML.Model)::Nothing
    _, formula, math = _parse_rule_formula(rule, model_SBML, libsbml_model; algebraic_rule=true)
    # As an algebraic rule gives the equation for a specie/parameter/compartment variable
    # is just a place-holder name
    if isempty(model_SBML.algebraic_rules)
        id = "1"
    else
        id = maximum(keys(model_SBML.algebraic_rules)) * "1"
    end
    have_rateOf = "rateOf" in math.fns
    model_SBML.algebraic_rules[id] = AlgebraicRuleSBML("0 ~ " * formula, math.math_idents, have_rateOf)
    return nothing
end

function _parse_rule_formula(rule::SBMLRule, model_SBML::ModelSBML, libsbml_model::SBML.Model; assignment_rule::Bool=false, rate_rule::Bool=false, algebraic_rule::Bool=false)::Tuple{String, String, MathSBML}
    @assert any([assignment_rule, algebraic_rule, rate_rule]) "No rule set to be parsed"

    variable = algebraic_rule ? "" : rule.variable
    rate_rule == true && push!(model_SBML.rate_rule_variables, variable)
    assignment_rule == true && push!(model_SBML.assignment_rule_variables, variable)

    math_expression = parse_math(rule.math, libsbml_model)

    formula = process_SBML_str_formula(math_expression.formula, model_SBML, libsbml_model; check_scaling=true, assignment_rule=assignment_rule, rate_rule=rate_rule, algebraic_rule=algebraic_rule, variable=variable)
    formula = isempty(formula) ? "0.0" : formula
    return variable, formula, math_expression
end

function _add_rule_info!(specie::SpecieSBML, formula::String, have_ridents::Bool, have_rateOf::Bool; assignment_rule::Bool=false, rate_rule::Bool=false)::Nothing
    specie.assignment_rule = assignment_rule
    specie.rate_rule = rate_rule
    if isempty(specie.initial_value) && !rate_rule
        specie.initial_value = formula
    end
    # If unit amount adjust for SBML equations per standard being given in conc.
    formula = _adjust_for_unit(formula, specie)
    specie.formula = formula
    specie.has_reaction_ids = have_ridents
    specie.has_rateOf = have_rateOf
    return nothing
end
function _add_rule_info!(variable::Union{ParameterSBML, CompartmentSBML}, formula::String, have_ridents::Bool, have_rateOf::Bool; assignment_rule::Bool=false, rate_rule::Bool=false)::Nothing
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
    variable.has_rateOf = have_rateOf
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

function identify_algebraic_rule_variables!(model_SBML::ModelSBML)::Nothing
    isempty(model_SBML.algebraic_rules) && return nothing
    for rule in values(model_SBML.algebraic_rules)
        # First species are checked, as due to edge cases several species can be valid,
        # in this case they must be compared between each other, see below
        algebraic_specie = _get_algebraic_species(rule.idents, model_SBML)
        if !isnothing(algebraic_specie)
            algebraic_specie.algebraic_rule = true
            push!(model_SBML.algebraic_rule_variables, algebraic_specie.name)
            continue
        end

        for ident in rule.idents
            variable = _get_model_variable(ident, model_SBML)
            variable isa SpecieSBML && continue
            _isassigned(variable) == true && continue
            variable.algebraic_rule = true
            push!(model_SBML.algebraic_rule_variables, variable.name)
            break
        end
    end
    return nothing
end

function _get_algebraic_species(idents, model_SBML)::Union{SpecieSBML, Nothing}
    valid_species = SpecieSBML[]
    for ident in idents
        variable = _get_model_variable(ident, model_SBML)
        !(variable isa SpecieSBML) && continue
        _isassigned(variable) == true && continue

        in_reactions = variable.name in model_SBML.species_in_reactions
        boundary_condition = variable.boundary_condition
        if !(in_reactions == false || boundary_condition)
            continue
        end
        push!(valid_species, variable)
    end
    length(valid_species) == 0 && return nothing
    length(valid_species) == 1 && return valid_species[1]
    # A rare edge case can happen where two species are in a sense valid. Here the rule
    # variable becomes the one not involved in a reaction, and if both are not involved
    # reactions, the one without a boundary condition is rule variable. Just SBML stuff...
    not_in_reactions = [s.name ∉ model_SBML.species_in_reactions for s in valid_species]
    if !all(not_in_reactions)
        return valid_species[findfirst(not_in_reactions)]
    end
    boundary_cond = [s.boundary_condition == false for s in valid_species]
    if !all(boundary_cond)
        return valid_species[findfirst(boundary_cond)]
    end
    throw(SBMLSupport("Cannot process algebraic rule, two species are equally valid"))
end

function _isassigned(variable::Union{SpecieSBML, ParameterSBML, CompartmentSBML})::Bool
    variable.constant && return true
    variable.rate_rule && return true
    variable.assignment_rule && return true
    return false
end
