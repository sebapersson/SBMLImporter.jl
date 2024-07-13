function _replace_variable(formula::AbstractString, to_replace::String,
                           replace_with::String)::String
    if !occursin(to_replace, formula)
        return formula
    end
    to_replace_r = Regex("(\\b" * to_replace * "\\b)")
    return replace(formula, to_replace_r => replace_with)
end

function _process_formula(math_expression::MathSBML, model_SBML::ModelSBML;
                          rate_rule::Bool = false, assignment_rule::Bool = false,
                          algebraic_rule::Bool = false, variable = "")::String
    @unpack formula, user_fns = math_expression
    formula = insert_functions(formula, model_SBML.functions, user_fns)
    formula = _replace_variable(formula, "time", "t")
    if _has_piecewise(formula)
        if algebraic_rule
            throw(SBMLSupport("Piecewise in algebraic rules is not supported"))
        end
        formula = piecewise_to_ifelse(formula)
        if (assignment_rule || rate_rule) && !isempty(variable)
            push!(model_SBML.variables_with_piecewise, variable)
        end
    end

    # Per the standard SBML formulas are given in concentration, therefore species with
    # unit amount must be rescaled.
    for (specie_id, specie) in model_SBML.species
        specie.unit == :Concentration && continue
        specie.only_substance_units && continue
        specie_id ∉ math_expression.math_idents && continue
        c = specie.compartment
        formula = _replace_variable(formula, specie_id, "(" * specie_id * "/" * c * ")")
    end
    return formula
end

function _has_piecewise(formula::String)::Bool
    return any(occursin.(["piecewise", "piecewise2", "piecewise4"] .* '(', formula))
end

function _get_model_as_str(path_SBML::String, model_as_string::Bool)::String
    if model_as_string == false
        if !isfile(path_SBML)
            throw(SBMLSupport("$path_SBML is not the path to a SBML file"))
        end
        f = open(path_SBML, "r")
        model_str = read(f, String)
        close(f)
        return model_str
    else
        return path_SBML
    end
end

function _has_time(formula::String)::Bool
    _formula = _replace_variable(formula, "t", "")
    return formula != _formula
end

function _contains_time_or_species(formula::String, model_SBML)::Bool
    _has_time(formula) == true && return true
    return any(occursin.(keys(model_SBML.species), formula))
end

function _trim_paranthesis(formula::String)::String
    formula = replace(formula, " " => "")
    length(formula) == 1 && return formula
    for _ in 1:100
        if (formula[1] == '(' && formula[end] == ')') &&
           (formula[2] == '(' && formula[end - 1] == ')')
            formula = formula[2:(end - 1)]
        else
            return formula
        end
    end
    return formula
end

function inline_assignment_rules!(model_SBML::ModelSBML)::Nothing
    for (reaction_id, reaction) in model_SBML.reactions
        if reaction.has_assignment_rule_variable == false
            continue
        end
        # Assignment rule can be nested so must go via recursion here, and
        # break when after one iteration kinetic-math formula has not changed
        while true
            _kinetic_math = reaction.kinetic_math

            for variable in model_SBML.assignment_rule_variables
                if !occursin(variable, reaction.kinetic_math)
                    continue
                end

                if haskey(model_SBML.species, variable)
                    formula = model_SBML.species[variable].formula
                elseif haskey(model_SBML.parameters, variable)
                    formula = model_SBML.parameters[variable].formula
                elseif haskey(model_SBML.compartments, variable)
                    formula = model_SBML.compartments[variable].formula
                end
                reaction.kinetic_math = _replace_variable(reaction.kinetic_math, variable,
                                                          formula)
            end

            if reaction.kinetic_math == _kinetic_math
                break
            end
        end

        # Need to inline into rate-rules
        for rule_id in model_SBML.rate_rule_variables
            if haskey(model_SBML.compartments, rule_id)
                raterule = model_SBML.compartments[rule_id]
            elseif haskey(model_SBML.parameters, rule_id)
                raterule = model_SBML.parameters[rule_id]
            elseif haskey(model_SBML.species, rule_id)
                raterule = model_SBML.species[rule_id]
            end
            # Assignment rule can be nested so must go via recursion here
            while true
                _kinetic_math = deepcopy(raterule.formula)
                for variable in model_SBML.assignment_rule_variables
                    if !occursin(variable, raterule.formula)
                        continue
                    end

                    if haskey(model_SBML.species, variable)
                        formula = model_SBML.species[variable].formula
                    elseif haskey(model_SBML.parameters, variable)
                        formula = model_SBML.parameters[variable].formula
                    elseif haskey(model_SBML.compartments, variable)
                        formula = model_SBML.compartments[variable].formula
                    end
                    raterule.formula = _replace_variable(reaction.kinetic_math, variable,
                                                         formula)
                end
                if raterule.formula == _kinetic_math
                    break
                end
            end
        end
    end

    # As assignment rule variables have been inlined they can be removed from the model
    filter!(isempty, model_SBML.assignment_rule_variables)

    return nothing
end

"""
    is_number(x::String)::Bool

    Check if a string x is a number (Float) taking scientific notation into account.
"""
function is_number(x::Union{AbstractString, SubString{String}})::Bool
    re1 = r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)$" # Picks up scientific notation
    re2 = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$"
    return (occursin(re1, x) || occursin(re2, x))
end

function _parse_bool(x::Union{Nothing, Bool})::Bool
    isnothing(x) && return false
    return x
end

function _parse_variable(x::Union{Nothing, String, Real}; default::String = "")::String
    isnothing(x) && return default
    return string(x)
end

function _get_times_appear(id::String, formula::AbstractString;
                           isfunction::Bool = true)::Integer
    pattern = isfunction ? Regex("\\b" * id * "\\(") : Regex("\\b" * id)
    return count(pattern, formula)
end

function _adjust_for_unit(formula::String, specie::SpecieSBML)::String
    if specie.unit == :Amount && specie.only_substance_units == false
        formula = '(' * formula * ") * " * specie.compartment
    end
    return formula
end

function _is_event_assigned(variable::Union{Nothing, String}, model_SBML::ModelSBML)::Bool
    isempty(model_SBML.events) && return false
    isnothing(variable) && return false
    for e in values(model_SBML.events)
        for formula in e.formulas
            formula_lhs = split(formula, "=")[1]
            occursin(variable, formula_lhs) && return true
        end
    end
    return false
end

function _get_specieref(variable::String,
                        libsbml_model::SBML.Model)::Union{Nothing, SBML.SpeciesReference}
    for reactions in values(libsbml_model.reactions)
        for reactant in reactions.reactants
            reactant.id == variable && return reactant
        end
        for product in reactions.products
            product.id == variable && return product
        end
    end
    return nothing
end

function _is_model_variable(variable::String, model_SBML::ModelSBML)::Bool
    haskey(model_SBML.species, variable) && return true
    haskey(model_SBML.parameters, variable) && return true
    haskey(model_SBML.compartments, variable) && return true
    return false
end

function _get_model_variable(variable::String, model_SBML::ModelSBML)
    haskey(model_SBML.species, variable) && return model_SBML.species[variable]
    haskey(model_SBML.parameters, variable) && return model_SBML.parameters[variable]
    haskey(model_SBML.compartments, variable) && return model_SBML.compartments[variable]
end

function _add_ident_info!(variable::VariableSBML, math_expression::MathSBML,
                          model_SBML::ModelSBML)::Nothing
    @unpack has_reaction_ids, has_rateOf, has_specieref = variable
    reaction_ids = math_expression.has_reaction_ids
    rateOf = math_expression.has_rateOf
    specieref = _has_specieref(math_expression.math_idents, model_SBML)
    variable.has_reaction_ids = _update_only_if_true(has_reaction_ids, reaction_ids)
    variable.has_rateOf = _update_only_if_true(has_rateOf, rateOf)
    variable.has_specieref = _update_only_if_true(has_specieref, specieref)
    return nothing
end

function _update_only_if_true(cond_old::Bool, cond_new::Bool)::Bool
    return cond_old == true ? true : cond_new
end

function _has_assignment_rule_ident(idents::Vector{String}, model_SBML::ModelSBML)::Bool
    for ident in idents
        ident in model_SBML.assignment_rule_variables && return true
    end
    return false
end

function _has_specieref(idents::Vector{String}, model_SBML::ModelSBML)::Bool
    for ident in idents
        ident in model_SBML.specie_reference_ids && return true
    end
    return false
end

function _adjust_assignment_rule_variables(formula::String, model_SBML::ModelSBML)::String
    for variables in [model_SBML.compartments, model_SBML.parameters]
        for (variable_id, variable) in variables
            if variable.assignment_rule == false
                continue
            end
            formula = _replace_variable(formula, variable_id, variable.formula)
        end
    end
    return formula
end

function _find_indices_outside_paranthesis(x::Char, formula::AbstractString;
                                           start_depth = 0)::Vector{Integer}
    out = Vector{Int64}(undef, 0)
    paranthesis_depth = start_depth
    for (i, char) in pairs(formula)
        if char == '('
            paranthesis_depth += 1
        end
        if char == ')'
            paranthesis_depth -= 1
        end
        if paranthesis_depth == 0 && char == x
            push!(out, i)
        end
    end
    return out
end

function _split_by_indicies(str::String, indices::Vector{<:Integer}; istart = 1,
                            iend = 0)::Vector{String}
    isempty(str[istart:(end - iend)]) && return String[]
    length(indices) == 0 && return [str[istart:(end - iend)]]

    out = Vector{String}(undef, length(indices) + 1)
    for (j, index) in pairs(indices)
        out[j] = str[istart:(index - 1)]
        istart = index + 1
    end
    out[end] = str[(indices[end] + 1):(end - iend)]
    return out
end

function _get_cf_scaling(specie::SpecieSBML, model_SBML::ModelSBML)::String
    if isempty(specie.conversion_factor) && isempty(model_SBML.conversion_factor)
        return ""
    elseif !isempty(specie.conversion_factor)
        return specie.conversion_factor
    else
        return model_SBML.conversion_factor
    end
end

function get_specie_reference_ids(libsbml_model::SBML.Model)::Vector{String}
    specie_reference_ids = String[]
    for r in values(libsbml_model.reactions)
        for reactant in r.reactants
            if isnothing(reactant.id)
                continue
            end
            push!(specie_reference_ids, reactant.id)
        end
        for product in r.products
            if isnothing(product.id)
                continue
            end
            push!(specie_reference_ids, product.id)
        end
    end
    return unique(specie_reference_ids)
end

function get_rule_variables!(model_SBML::ModelSBML)::Nothing
    rule_variables = Iterators.flatten((model_SBML.rate_rule_variables,
                                        model_SBML.assignment_rule_variables,
                                        model_SBML.algebraic_rule_variables))
    rule_variables = unique(rule_variables)
    filter!(x -> !haskey(model_SBML.generated_ids, x), rule_variables)
    for var in rule_variables
        push!(model_SBML.rule_variables, var)
    end
    # For downstrem processing, further filter generated ids from other rule list
    filter!(x -> !haskey(model_SBML.generated_ids, x), rule_variables)
    filter!(x -> !haskey(model_SBML.generated_ids, x), model_SBML.assignment_rule_variables)
    filter!(x -> !haskey(model_SBML.generated_ids, x), model_SBML.rate_rule_variables)
    return nothing
end
