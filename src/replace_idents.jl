function replace_ident!(model_SBML::ModelSBML, libsbml_model::SBML.Model, ident_replace::Symbol)::Nothing
    @assert ident_replace in [:rateOf, :specieref, :reactionid] "Unknown ident $ident_replace to replace"
    if ident_replace == :rateOf
        f_replace = _replace_rateOf
    elseif ident_replace == :specieref
        f_replace = _replace_specieref
    elseif ident_replace == :reactionid
        f_replace = _replace_reactionid
    end

    @unpack species, parameters, compartments = model_SBML
    for (_, variable) in Iterators.flatten((species, parameters, compartments))
        _has_ident(variable, ident_replace) == false && continue
        variable.formula = f_replace(variable.formula, model_SBML, libsbml_model)
        variable.initial_value = f_replace(variable.initial_value, model_SBML, libsbml_model)
    end
    for rule in values(model_SBML.algebraic_rules)
        ident_replace in [:specieref, :reactionid] && continue
        rule.has_rateOf == false && continue
        rule.formula = f_replace(rule.formula, model_SBML, libsbml_model)
    end
    for reaction in values(model_SBML.reactions)
        _has_ident(reaction, ident_replace) == false && continue
        reaction.kinetic_math = f_replace(reaction.kinetic_math, model_SBML, libsbml_model)
    end
    for event in values(model_SBML.events)
        if _has_ident(event, true, ident_replace)
            event.trigger = f_replace(event.trigger, model_SBML, libsbml_model)
        end
        _has_ident(event, false, ident_replace) == false && continue
        for (i, formula) in pairs(event.formulas)
            event.formulas[i] = f_replace(formula, model_SBML, libsbml_model)
        end
    end
    return nothing
end

function _replace_rateOf(formula::String, model_SBML::ModelSBML, ::SBML.Model)::String
    nrateOf_occurrences = _get_times_appear("rateOf", formula)
    for _ in 1:nrateOf_occurrences
        function_call = _extract_function_call("rateOf", formula)
        rateOf_insert = _get_rateOf_insert(function_call, model_SBML)
        formula = replace(formula, function_call => rateOf_insert; count=1)
    end
    return formula
end

function _replace_specieref(formula::String, model_SBML::ModelSBML, libsbml_model::SBML.Model)::String
    # Sometimes specie_reference_ids appear in rules, kinetic-math..., where the reference
    # id it not assigned by any rules. In this case the reference id should be replaced by
    # its stoichemetry, or if provided via an initial assignment, its initial assignment
    # value
    for specie_reference_id in model_SBML.specie_reference_ids
        isnothing(specie_reference_id) && continue
        specie_reference_id in model_SBML.libsbml_rule_variables && continue

        if !haskey(libsbml_model.initial_assignments, specie_reference_id)
            haskey(libsbml_model.species, specie_reference_id) && continue
            specie_reference = _get_specieref(specie_reference_id, libsbml_model)
            S = string(specie_reference.stoichiometry)
            formula = _replace_variable(formula, specie_reference.id, S)
            continue
        end
        # Initial assignments might not map to any model variables, so reparsing of
        # formula is needed.
        assignment = libsbml_model.initial_assignments[specie_reference_id]
        s0 = _parse_assignment_formula(specie_reference_id, assignment, model_SBML, libsbml_model)
        s0 = "(" * s0.formula * ")"
        formula = _replace_variable(formula, specie_reference_id, s0)
    end
    return formula
end

function _replace_reactionid(formula::String, model_SBML::ModelSBML, ::SBML.Model)::String
    if !any(occursin.(keys(model_SBML.reactions), formula))
        return formula
    end
    for (reaction_id, reaction) in model_SBML.reactions
        formula = _replace_variable(formula, reaction_id, reaction.kinetic_math)
    end
    return formula
end

function _get_rateOf_insert(function_call::String, model_SBML::ModelSBML)::String
    arg = _extract_args_insert(function_call, true)[1]

    is_number(arg) && return "0.0"
    if haskey(model_SBML.parameters, arg)
        p = model_SBML.parameters[arg]
        p.constant == true && return "0.0"
        p.rate_rule == true && return p.formula
    end
    # If we get here rateOf correspond to a specie, but, the specie might be scaled in
    # earlier steps with compartment. If this is the case, specie is extracted first, and
    # corresponding kinetic_math is inserted.
    # If specie unit is amount, as SBML formulas are given in conc., a scaling is performed
    specie = _get_specie_rateOf(arg, model_SBML)
    if specie.unit == :Amount && specie.only_substance_units == false
        return "(" * specie.formula * ") / " * specie.compartment
    else
        return specie.formula
    end
end

function _get_specie_rateOf(arg::String, model_SBML::ModelSBML)::SpecieSBML
    if haskey(model_SBML.species, arg)
        return model_SBML.species[arg]
    end
    arg = filter(x -> x ∉ ['(', ')'], arg)
    arg = occursin('/', arg) ? arg[1:(findfirst(x -> x == '/', arg) - 1)] : arg
    arg = occursin('*', arg) ? arg[1:(findfirst(x -> x == '*', arg) - 1)] : arg
    @assert haskey(model_SBML.species, arg) "rateOf arg $arg is not a model specie"
    return model_SBML.species[arg]
end

function _has_ident(variable::Union{VariableSBML, ReactionSBML}, ident_replace::Symbol)::Bool
    if ident_replace == :rateOf
        return variable.has_rateOf
    elseif ident_replace == :specieref
        return variable.has_specieref
    else
        return variable.has_reaction_ids
    end
end
function _has_ident(event::EventSBML, trigger::Bool, ident_replace::Symbol)::Bool
    if ident_replace == :rateOf
        trigger == true && return event.has_rateOf_trigger
        trigger == false && return event.has_rateOf_assignments
    elseif ident_replace == :specieref
        trigger == true && return event.has_specieref_trigger
        trigger == false && return event.has_specieref_assignments
    else
        trigger == true && return event.has_reaction_ids_trigger
        trigger == false && return event.has_reaction_ids_assignments
    end
end
