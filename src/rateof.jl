function replace_rateOf!(model_SBML::ModelSBML)::Nothing
    variables = Iterators.flatten((model_SBML.species, model_SBML.parameters))
    for (_, variable) in variables
        variable.has_rateOf == false && continue
        variable.formula = _replace_rateOf(variable.formula, model_SBML)
        variable.initial_value = _replace_rateOf(variable.initial_value, model_SBML)
    end
    for rule in values(model_SBML.algebraic_rules)
        rule.has_rateOf == false && continue
        rule.formula = _replace_rateOf(rule.formula, model_SBML)
    end
    for reaction in values(model_SBML.reactions)
        reaction.has_rateOf == false && continue
        reaction.kinetic_math = _replace_rateOf(reaction.kinetic_math, model_SBML)
    end
    for event in values(model_SBML.events)
        if event.has_rateOf_trigger
            event.trigger = _replace_rateOf(event.trigger, model_SBML)
        end
        event.has_rateOf_assignments == false && continue
        for (i, formula) in pairs(event.formulas)
            event.formulas[i] = _replace_rateOf(formula, model_SBML)
        end
    end
    return nothing
end

function _replace_rateOf(formula::String, model_SBML::ModelSBML)::String
    nrateOf_occurrences = _get_times_appear("rateOf", formula)
    for _ in 1:nrateOf_occurrences
        function_call = _extract_function_call("rateOf", formula)
        rateOf_insert = _get_rateOf_insert(function_call, model_SBML)
        formula = replace(formula, function_call => rateOf_insert; count=1)
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
