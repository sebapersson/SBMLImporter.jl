# Handles piecewise functions that are to be redefined with ifelse speciements in the model
# equations to allow MKT symbolic calculations.
function piecewise_to_ifelse(formula::String)::String
    return insert_functions(formula, PIECEWISE_FN, PIECEWISE_FN_NAMES)
end

function time_dependent_ifelse_to_bool!(model_SBML::ModelSBML)::Nothing
    variables = Iterators.flatten((model_SBML.species, model_SBML.parameters))
    for (_, variable) in variables
        !occursin("ifelse", variable.formula) && continue
        variable.formula = _time_dependent_ifelse_to_bool(variable.formula, model_SBML)
    end
    return nothing
end

function _time_dependent_ifelse_to_bool(formula::String, model_SBML::ModelSBML;
                                        ifelse_parameter_names = String[])::String
    if !occursin("ifelse", formula)
        return formula
    end
    ifelse_with_time::Bool = false
    nifelse = _get_times_appear("ifelse", formula)
    ifelse_calls = _extract_function_calls("ifelse", formula)
    @assert length(ifelse_calls)==nifelse "Error in ifelse to bool parsing"

    for ifelse_call in ifelse_calls
        # In some cases the same ifelse-call can appear in several equations
        if haskey(model_SBML.ifelse_bool_expressions, ifelse_call)
            ifelse_with_time = true
            formula_bool = model_SBML.ifelse_bool_expressions[ifelse_call]
            formula = replace(formula, ifelse_call => formula_bool)
            break
        end

        # If !=, ==, false, true, ifelse is in condition rewriting to event is not possible.
        # This rarely happens outside SBML test-suite
        condition, arg1, arg2 = _extract_args_insert(ifelse_call)
        condition = _trim_paranthesis(condition)
        if any(occursin.(["!=", "==", "false", "true", "ifelse"], condition))
            continue
        end

        lhs_condition, rhs_condition, operator = _split_condition(condition)
        time_in_rhs = _has_time(rhs_condition)
        time_in_lhs = _has_time(lhs_condition)
        time_in_lhs == false && time_in_rhs == false && continue
        @assert time_in_rhs!=time_in_lhs "Error with time in both condition sides"
        side_activated = _get_side_activated_with_time(lhs_condition, rhs_condition,
                                                       operator, time_in_rhs)

        bool_name = _get_name_bool_piecewise(ifelse_parameter_names)
        formula_bool = _template_bool_picewise(bool_name, arg1, arg2, side_activated)
        formula = replace(formula, ifelse_call => formula_bool; count = 1)

        model_SBML.ifelse_bool_expressions[ifelse_call] = formula_bool
        model_SBML.parameters[bool_name] = ParameterSBML(bool_name, true, "0.0", "", false,
                                                         false, false, false, false, false)
        model_SBML.events[bool_name] = _ifelse_to_event(bool_name, condition,
                                                        side_activated)

        ifelse_with_time = true
        break
    end

    if ifelse_with_time == true
        formula = _time_dependent_ifelse_to_bool(formula, model_SBML;
                                                 ifelse_parameter_names = ifelse_parameter_names)
    end
    return formula
end

function _get_side_activated_with_time(lhs_condition::String, rhs_condition::String,
                                       operator::String, time_in_rhs::Bool)::String
    sign_time = time_in_rhs ? _get_sign_time(rhs_condition) : _get_sign_time(lhs_condition)
    # Example : if we have -t > -1 then sign_time = -1, and when time increases from
    # t0=0 we have with time a transition from true -> false, which means that in the
    # ifelse the right side is activated with time
    if operator in ["<", "≤", "<="]
        side_if_rhs = sign_time == 1 ? "left" : "right"
    end
    if operator in [">", "≥", ">="]
        side_if_rhs = sign_time == 1 ? "right" : "left"
    end
    side_if_lhs = side_if_rhs == "left" ? "right" : "left"
    if time_in_rhs
        return side_if_rhs
    else
        return side_if_lhs
    end
end

function _split_condition(formula::String)::Tuple{String, String, String}
    gt_applys = [">", "≥", ">="]
    lt_applys = ["<", "≤", "<="]
    igt = findfirst(x -> occursin(x, formula), gt_applys)
    ilt = findfirst(x -> occursin(x, formula), lt_applys)
    @assert !all(isnothing.([igt, ilt])) "Error splitting ifelse condition"
    operator = isnothing(igt) ? lt_applys[ilt] : gt_applys[igt]
    lhs, rhs = string.(split(formula, operator))
    return lhs, rhs, operator
end

function _get_name_bool_piecewise(ifelse_parameter_names::Vector{String})::String
    j = 1
    while true
        parameter_name = "__parameter_ifelse" * string(j)
        !(parameter_name in ifelse_parameter_names) && break
        j += 1
    end
    parameter_name = "__parameter_ifelse" * string(j)
    push!(ifelse_parameter_names, parameter_name)
    return parameter_name
end

function _template_bool_picewise(bool_name::String, ifelse_arg1::String,
                                 ifelse_arg2::String, side_activated::String)::String
    activated = side_activated == "left" ? ifelse_arg1 : ifelse_arg2
    deactivated = side_activated == "left" ? ifelse_arg2 : ifelse_arg1
    formula = "((1 - 1" * bool_name * ") * (" * deactivated * ") + " *
              bool_name * "*(" * activated * "))"
    return formula
end

function _get_sign_time(formula::String)::Int64
    formula = replace(formula, " " => "")
    formula = _trim_paranthesis(formula)
    _formula = _find_term_with_t(formula)
    @assert !isempty(_formula) "In $formula for condition in piecewise cannot identify which term time appears in."

    _formula = replace(_formula, "(" => "", ")" => "")
    if _formula == "t"
        return 1
    elseif length(_formula) ≥ 2 && _formula[1:2] == "-t"
        return -1
    elseif length(_formula) ≥ 3 && _formula[1:3] == "--t"
        return 1
    elseif length(_formula) ≥ 3 && (_formula[1:3] == "+-t" || _formula[1:3] == "-+t")
        return -1
    end

    # If '-' appears anywhere after the cases above we might be able to infer direction,
    # but infering in this situation is hard! - so throw an error as the user should
    # be able to write condition in a more easy manner (avoid several sign changing minus signs)
    !occursin('-', _formula) && return 1
    str_write = "For piecewise with time in condition we cannot infer direction for $formula, that is if the condition value increases or decreases with time. This happens if the formula contains a minus sign in the term where t appears."
    throw(SBMLSupport(str_write))
end

function _find_term_with_t(formula::String)::String
    formula == "t" && return formula
    istart, parenthesis_level, term = 1, 0, ""
    for (i, char) in pairs(formula)
        if i == 1 && char in ['+', '-']
            continue
        end
        parenthesis_level = char == '(' ? parenthesis_level + 1 : parenthesis_level
        parenthesis_level = char == ')' ? parenthesis_level - 1 : parenthesis_level
        if !(parenthesis_level == 0 && (char in ['+', '-'] || i == length(formula)))
            continue
        end
        if formula[i - 1] in ['+', '-'] && i != length(formula)
            continue
        end
        term = formula[istart:i]
        if _has_time(term) == true || (i == length(formula) && char == 't')
            return term[end] in ['+', '-'] ? term[1:(end - 1)] : term
        end
        istart = i
    end
    return ""
end

function _ifelse_to_event(id::String, condition::String, side_activated)::EventSBML
    assignments = [id * " = 1.0"]
    # When triggered we change the bool variable from 0 to 1. If side_activated = "right"
    # we activate the event when the ifelse condition goes from true to false. Reverse for
    # side_activated = "left". Thus, for side_activated = "right" we need to invert the
    # condition, as otherwise we mess up callback initialisation where the bool variable
    # is set to 1 if the condition is true
    if side_activated == "right"
        if any(occursin.(["<", "≤", "<="], condition))
            condition = replace(condition, r"≤|<=|<" => "≥")
        else
            condition = replace(condition, r"≥|>=|>" => "≤")
        end
    end
    event = EventSBML(id, condition, assignments, false, false, false, false, false, false,
                      false, true)
    return event
end
