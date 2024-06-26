function parse_initial_assignments!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (assign_id, assignment) in libsbml_model.initial_assignments
        math_expression = _parse_assignment_formula(assign_id, assignment, model_SBML, libsbml_model)

        if _is_model_variable(assign_id, model_SBML) == true
            _add_assignment_info!(model_SBML, assign_id, math_expression)
            continue
        end
        @warn "Initial assignment creates new specie $(assign_id). Happens when " *
              "$(assign_id) does not correspond to any model specie, parameter, or compartment."
        model_SBML.species[assign_id] = SpecieSBML(assign_id, false, false, math_expression.formula, "", "1.0", "", :Amount, false, false, false, false, false, false, false)
        _add_ident_info!(model_SBML.species[assign_id], math_expression, model_SBML)
    end

    # The initial assignment can depend on other species, but this cannot be handled by
    # Catalyst.jl. Therefore, recursively unnest species.
    for specie_id in keys(model_SBML.species)
        !haskey(libsbml_model.initial_assignments, specie_id) && continue
        _unnest_initial_assignment!(model_SBML, specie_id)
    end
    # If a specie formula depend on assignment rule computed compartments, these must
    # be inlined as in the initial-value mapping assignment rules are not identified.
    # Holds for all species
    for specie in values(model_SBML.species)
        specie.initial_value = _adjust_assignment_rule_variables(specie.initial_value, model_SBML)
    end
    return nothing
end

function _parse_assignment_formula(assignment_id::String, assignment, model_SBML::ModelSBML, libsbml_model::SBML.Model)::MathSBML
    # Initial assignment applies at t = 0.0
    math_expression = parse_math(assignment, libsbml_model)
    formula = _process_formula(math_expression, model_SBML; variable=assignment_id)
    formula = _replace_variable(formula, "t", "0.0")
    math_expression.formula = formula
    return math_expression
end

function _add_assignment_info!(model_SBML::ModelSBML, assign_id::String, math_expression::MathSBML)::Nothing
    variable = _get_model_variable(assign_id, model_SBML)
    _add_ident_info!(variable, math_expression, model_SBML)
    formula = math_expression.formula
    if haskey(model_SBML.species, variable.name)
        variable.initial_value = _adjust_for_unit(formula, variable)
        return nothing
    end

    if assign_id in model_SBML.rate_rule_variables
        variable.initial_value = formula
        return nothing
    end
    # The variable can now be non-constant (depend on species). This can be handled by
    # treating it as an assignment rule
    variable.formula = formula
    variable.initial_value = formula
    variable.assignment_rule = true
    push!(model_SBML.assignment_rule_variables, assign_id)
    return nothing
end

function _unnest_initial_assignment!(model_SBML::ModelSBML, assign_id::String)::Nothing
    # rateOf expression are handled later by the importer (everything must be parsed)
    variable = _get_model_variable(assign_id, model_SBML)
    if variable.has_rateOf == true
        return nothing
    end
    formula = variable.initial_value

    for i in 1:100
        _formula = deepcopy(formula)
        for (specie_id, specie) in model_SBML.species
            if specie_id == assign_id
                continue
            end
            formula = _replace_variable(formula, specie_id, specie.initial_value)
        end
        if formula == _formula
            break
        end
        @assert i != 100 "Got stuck in recursion while unnesting initial assignment"
    end

    variable.initial_value = formula
    return nothing
end
