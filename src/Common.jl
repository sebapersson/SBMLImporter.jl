"""
    replace_variable(formula::T, to_replace::String, replace_with::String)::T where T<:AbstractString

In formula, replaces to_replace with replace_with. Exact match is required, so if to_replace=time1
and replace_with=1 while formula = time * 2 nothing is replaced.
"""
function replace_variable(formula::T, to_replace::String,
                          replace_with::String)::T where {T <: AbstractString}
    if !occursin(to_replace, formula)
        return formula
    end

    _to_replace = Regex("(\\b" * to_replace * "\\b)")
    return replace(formula, _to_replace => replace_with)
end

"""
    process_SBML_str_formula(formula::T, model_SBML::ModelSBML, libsbml_model::SBML.Model;
                             check_scaling::Bool=false, rate_rule::Bool=false)::T where T<:AbstractString

Processes a string formula by inserting SBML functions, rewriting piecewise to ifelse, and scaling species.
"""
function process_SBML_str_formula(formula::T, model_SBML::ModelSBML,
                                  libsbml_model::SBML.Model;
                                  check_scaling::Bool = false,
                                  rate_rule::Bool = false,
                                  assignment_rule::Bool = false,
                                  algebraic_rule::Bool = false,
                                  scale_rateof::Bool = true,
                                  initial_assignment::Bool = false,
                                  variable = "")::T where {T <: AbstractString}
    _formula = insert_functions(formula, model_SBML.functions)
    if algebraic_rule && tmp_has_piecewise(formula)
        throw(SBMLSupport("Piecewise in algebraic rules is not supported"))
    end
    if initial_assignment && haskey(model_SBML.species, variable) && tmp_has_piecewise(formula)
        throw(SBMLSupport("Piecewise in initial assignment is not supported"))
    end
    if tmp_has_piecewise(_formula)
        _formula = piecewise_to_ifelse(_formula)
        if assignment_rule || rate_rule
            !isempty(variable) && push!(model_SBML.variables_with_piecewise, variable)
        end
    end
    _formula = replace_variable(_formula, "time", "t") # Sometimes t is decoded as time

    # SBML equations are given in concentration, in case an amount specie appears in the equation scale with the
    # compartment in the formula every time the species appear
    for (specie_id, specie) in model_SBML.species
        if !occursin(specie_id, _formula)
            continue
        end
        if check_scaling == false
            continue
        end
        if specie.unit == :Concentration || specie.only_substance_units == true
            continue
        end
        # TODO: Remember for rateOf parsing
        if scale_rateof == false && occursin("rateOf", _formula)
            continue
        end

        compartment = specie.compartment
        _formula = replace_variable(_formula, specie_id,
                                    "(" * specie_id * "/" * compartment * ")")
    end

    # Replace potential expressions given in initial assignment and that appear in stoichemetric experssions
    # of reactions (these are not species, only math expressions that should be replaced)
    for id in keys(libsbml_model.initial_assignments)
        # In case ID does not occur in stoichemetric expressions
        if isempty(libsbml_model.reactions)
            continue
        end
        if id ∉ model_SBML.specie_reference_ids
            continue
        end
        if !haskey(model_SBML.species, id) && rate_rule == false
            continue
        end
        if isnothing(id)
            continue
        end
        # Do not rewrite is stoichemetric is controlled via event
        if !isempty(model_SBML.events) &&
           any(occursin.(id, reduce(vcat, [e.formulas for e in values(model_SBML.events)])))
            continue
        end
        if rate_rule == false
            _formula = replace_variable(_formula, id,
                                        "(" * model_SBML.species[id].initial_value * ")")
        else
            math_expression = parse_math(libsbml_model.initial_assignments[id], libsbml_model)
            _formula = replace_variable(_formula, id, "(" * math_expression.formula * ")")
        end
    end

    # Sometimes we have a stoichemetric expression appearing in for example rule expressions, etc..., but it does not
    # have any initial assignment, or rule assignment. In this case the reference should be replaced with its
    # corresponding stoichemetry
    if any(occursin.(model_SBML.specie_reference_ids, _formula))
        for (_, reaction) in libsbml_model.reactions
            specie_references = Iterators.flatten(([reactant
                                                    for reactant in reaction.reactants],
                                                   [product
                                                    for product in reaction.products]))
            for specie_reference in specie_references
                if isnothing(specie_reference.id)
                    continue
                end
                if haskey(libsbml_model.species, specie_reference.id)
                    continue
                end
                if haskey(libsbml_model.initial_assignments, specie_reference.id)
                    continue
                end
                if specie_reference.id ∈ [rule isa SBML.AlgebraicRule ? "" : rule.variable
                    for rule in libsbml_model.rules]
                    continue
                end
                _formula = replace_variable(_formula, specie_reference.id,
                                            string(specie_reference.stoichiometry))
            end
        end
    end

    return _formula
end

# TODO: Not needed when start tracking fn
function tmp_has_piecewise(formula)::Bool
    occursin("piecewise(", formula) && return true
    occursin("piecewise2(", formula) && return true
    occursin("piecewise4(", formula) && return true
    return false
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
    _formula = replace_variable(formula, "t", "")
    return formula != _formula
end

function _trim_paranthesis(formula::String)::String
    length(formula) == 1 && return formula
    for _ in 1:100
        if (formula[1] == '(' && formula[end] == ')') &&
           (formula[2] == '(' && formula[end-1] == ')')
            formula = formula[2:end-1]
        else
            return formula
        end
    end
    return formula
end

function replace_rateOf!(model_SBML::ModelSBML)::Nothing
    for (parameter_id, parameter) in model_SBML.parameters
        parameter.formula = replace_rateOf(parameter.formula, model_SBML)
        parameter.initial_value = replace_rateOf(parameter.initial_value, model_SBML)
    end
    for (specie_id, specie) in model_SBML.species
        specie.formula = replace_rateOf(specie.formula, model_SBML)
        specie.initial_value = replace_rateOf(specie.initial_value, model_SBML)
    end
    for (event_id, event) in model_SBML.events
        for (i, formula) in pairs(event.formulas)
            event.formulas[i] = replace_rateOf(formula, model_SBML)
        end
        event.trigger = replace_rateOf(event.trigger, model_SBML)
    end
    for (rule_id, rule) in model_SBML.algebraic_rules
        rule.formula = replace_rateOf(rule.formula, model_SBML)
    end
    for (reaction_id, reaction) in model_SBML.reactions
        reaction.kinetic_math = replace_rateOf(reaction.kinetic_math, model_SBML)
    end

    return nothing
end

function replace_rateOf(_formula::T,
                        model_SBML::ModelSBML)::String where {
                                                              T <: Union{<:AbstractString,
                                                                    <:Real}}
    formula = string(_formula)
    if !occursin("rateOf", formula)
        return formula
    end

    # Invalid character problems
    formula = replace(formula, "≤" => "<=")
    formula = replace(formula, "≥" => ">=")

    # Find rateof expressions
    start_rateof = findall(i -> formula[i:(i + 6)] == "rateOf(", 1:(length(formula) - 6))
    end_rateof = [findfirst(x -> x == ')', formula[start:end]) + start - 1
                  for start in start_rateof]
    # Compenstate for nested paranthesis
    for i in eachindex(end_rateof)
        if any(occursin.(['*', '/'], formula[start_rateof[i]:end_rateof[i]]))
            end_rateof[i] += 1
        end
    end
    args = [formula[(start_rateof[i] + 7):(end_rateof[i] - 1)]
            for i in eachindex(start_rateof)]

    replace_with = Vector{String}(undef, length(args))
    for (i, arg) in pairs(args)

        # A constant parameter does not have a rate
        if haskey(model_SBML.parameters, arg) && model_SBML.parameters[arg].constant == true
            replace_with[i] = "0.0"
            continue
        end
        # A parameter via a rate-rule has a rate
        if haskey(model_SBML.parameters, arg) &&
           model_SBML.parameters[arg].rate_rule == true
            replace_with[i] = model_SBML.parameters[arg].formula
            continue
        end

        # A number does not have a rate
        if is_number(arg)
            replace_with[i] = "0.0"
            continue
        end

        # If specie is via a rate-rule we do not scale the state in the expression
        if haskey(model_SBML.species, arg) && model_SBML.species[arg].rate_rule == true
            replace_with[i] = model_SBML.species[arg].formula
            continue
        end

        # Default case, use formula for given specie, and if specie is given in amount
        # Here it might happen that arg is scaled with compartment, e.g. S / C thus
        # first the specie is extracted
        arg = filter(x -> x ∉ ['(', ')'], arg)
        arg = occursin('/', arg) ? arg[1:(findfirst(x -> x == '/', arg) - 1)] : arg
        arg = occursin('*', arg) ? arg[1:(findfirst(x -> x == '*', arg) - 1)] : arg
        specie = model_SBML.species[arg]
        scale_with_compartment = specie.unit == :Amount &&
                                 specie.only_substance_units == false
        if scale_with_compartment == true
            replace_with[i] = "(" * specie.formula * ") / " * specie.compartment
        else
            replace_with[i] = specie.formula
        end
    end

    formula_cp = deepcopy(formula)
    for i in eachindex(replace_with)
        formula = replace(formula,
                          formula_cp[start_rateof[i]:end_rateof[i]] => replace_with[i])
    end

    formula = replace(formula, "<=" => "≤")
    formula = replace(formula, ">=" => "≥")

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
                reaction.kinetic_math = replace_variable(reaction.kinetic_math, variable,
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
                    raterule.formula = replace_variable(reaction.kinetic_math, variable,
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

function _parse_variable(x::Union{Nothing, String, Real}; default::String="")::String
    isnothing(x) && return default
    return string(x)
end

function _get_times_appear(id::String, formula::AbstractString; isfunction::Bool=true)::Integer
    pattern = isfunction ? Regex("\\b" * id * "\\(") : Regex("\\b" * id)
    return count(pattern, formula)
end

function _adjust_for_unit(formula::String, specie::SpecieSBML)::String
    if specie.unit == :Amount && specie.only_substance_units == false
        formula = '(' * formula * ") * " * specie.compartment
    end
    return formula
end

# TODO: What is assigned to should be its own field
function _is_event_variable(variable::Union{Nothing, String}, model_SBML::ModelSBML)::Bool
    isnothing(variable) && return false
    for e in values(model_SBML.events)
        for formula in e.formulas
            occursin(variable, formula) && return true
        end
    end
    return false
end

function _which_specieref(variable::String, libsbml_model::SBML.Model)::Union{Nothing, SBML.SpeciesReference}
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
    @assert _is_model_variable(variable, model_SBML) "$variable not a SBML id"
    haskey(model_SBML.species, variable) && return model_SBML.species[variable]
    haskey(model_SBML.parameters, variable) && return model_SBML.parameters[variable]
    haskey(model_SBML.compartments, variable) && return model_SBML.compartments[variable]
end

function _has_assignment_rule_ident(idents::Vector{String}, model_SBML::ModelSBML)::Bool
    for ident in idents
        ident in model_SBML.assignment_rule_variables && return true
    end
    return false
end

function _has_reactionid_ident(idents::Vector{String}, libsbml_model::SBML.Model)::Bool
    for ident in idents
        haskey(libsbml_model.reactions, ident) && return true
    end
    return false
end

function _adjust_assignment_rule_variables(formula::String, model_SBML::ModelSBML; initial_assignment::Bool=false)::String
    for variables in [model_SBML.compartments, model_SBML.parameters]
        for (variable_id, variable) in variables
            if variable.assignment_rule == false
                continue
            end
            formula = replace_variable(formula, variable_id, variable.formula)
        end
    end
    return formula
end

function _find_indices_outside_paranthesis(x::Char, formula::AbstractString; start_depth=0)::Vector{Integer}
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

function _split_by_indicies(str::String, indices::Vector{<:Integer}; istart=1, iend=0)::Vector{String}
    isempty(str[istart:end-iend]) && return String[]
    length(indices) == 0 && return [str[istart:end-iend]]

    out = Vector{String}(undef, length(indices)+1)
    for (j, index) in pairs(indices)
        out[j] = str[istart:index-1]
        istart = index + 1
    end
    out[end] = str[(indices[end]+1):end-iend]
    return out
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
    return specie_reference_ids
end

function get_rule_variables!(model_SBML::ModelSBML)::Nothing
    _rule_variables = unique(vcat(model_SBML.rate_rule_variables,
                                  model_SBML.assignment_rule_variables,
                                  model_SBML.algebraic_rule_variables))
    filter!(x -> x ∉ keys(model_SBML.generated_ids), _rule_variables)

    for var in _rule_variables
        push!(model_SBML.rule_variables, var)
    end

    return nothing
end
