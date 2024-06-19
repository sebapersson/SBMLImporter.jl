function parse_initial_assignments!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (assign_id, assignment) in libsbml_model.initial_assignments
        formula, have_ridents = _parse_assignment_formula(assign_id, assignment, model_SBML, libsbml_model)

        if _is_model_variable(assign_id, model_SBML) == true
            _add_assignment_info!(model_SBML, assign_id, formula, have_ridents)
            continue
        end
        @warn "Initial assignment creates new specie $(assign_id). Happens when " *
              "$(assign_id) does not correspond to any model specie, parameter, or compartment."
        model_SBML.species[assign_id] = SpecieSBML(assign_id, false, false, formula, "", "1.0", "", :Amount, false, false, false, false, have_ridents)
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

function _parse_assignment_formula(assignment_id::String, assignment, model_SBML::ModelSBML, libsbml_model::SBML.Model)::Tuple{String, Bool}
    # Initial assignment applies at t = 0.0
    math_expression = parse_math(assignment, libsbml_model)
    have_ridents = _has_reactionid_ident(math_expression.math_idents, libsbml_model)
    formula = process_SBML_str_formula(math_expression.formula, model_SBML, libsbml_model; check_scaling=true, variable=assignment_id, initial_assignment=true)
    formula = replace_variable(formula, "t", "0.0")
    # TODO: This information lives in math-expression, and allows for early exit
    if any(occursin.(["and", "or"], formula))
        throw(SBMLSupport("and, or not supported in initial assignments"))
    end
    return formula, have_ridents
end

function _add_assignment_info!(model_SBML::ModelSBML, assign_id::String, formula::String, have_ridents::Bool)::Nothing
    variable = _get_model_variable(assign_id, model_SBML)
    variable.has_reaction_ids = have_ridents
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
    formula = variable.initial_value
    if occursin("rateOf", formula)
        return nothing
    end

    for i in 1:100
        _formula = deepcopy(formula)
        for (specie_id, specie) in model_SBML.species
            if specie_id == assign_id
                continue
            end
            formula = replace_variable(formula, specie_id, specie.initial_value)
        end
        if formula == _formula
            break
        end
        @assert i != 100 "Got stuck in recursion while unnesting initial assignment"
    end

    variable.initial_value = formula
    return nothing
end
