# Handles piecewise functions that are to be redefined with ifelse speciements in the model
# equations to allow MKT symbolic calculations.
function piecewise_to_ifelse(formula::String)::String
    return insert_functions(formula, PIECEWISE_FN, PIECEWISE_FN_NAMES)
end

function time_dependent_ifelse_to_bool!(model_SBML::ModelSBML)::Nothing
    variables = Iterators.flatten((model_SBML.species, model_SBML.parameters))
    ifelse_parameter_names = String[]
    for (_, variable) in variables
        !occursin("ifelse", variable.formula) && continue
        variable.formula = _time_dependent_ifelse_to_bool(variable.formula, model_SBML, ifelse_parameter_names)
    end
    return nothing
end

function _time_dependent_ifelse_to_bool(
        formula::String, model_SBML::ModelSBML, ifelse_parameter_names::Vector{String}
    )::String
    if !occursin("ifelse", formula)
        return formula
    end
    ifelse_with_time::Bool = false
    nifelse = _get_times_appear("ifelse", formula)
    ifelse_calls = _extract_function_calls("ifelse", formula)
    @assert length(ifelse_calls) == nifelse "Error in ifelse to bool parsing"

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
        @assert time_in_rhs != time_in_lhs "Error with time in both condition sides"
        side_activated = _get_side_activated_with_time(lhs_condition, rhs_condition, operator, time_in_rhs)

        bool_name = _get_name_bool_piecewise(ifelse_parameter_names)
        formula_bool = _template_bool_picewise(bool_name, arg1, arg2, side_activated)
        formula = replace(formula, ifelse_call => formula_bool; count = 1)

        model_SBML.ifelse_bool_expressions[ifelse_call] = formula_bool
        model_SBML.parameters[bool_name] = ParameterSBML(
            bool_name, true, "0.0", "", false, false, false, false, false, false
        )
        model_SBML.events[bool_name] = _ifelse_to_event(bool_name, condition, side_activated)

        ifelse_with_time = true
        break
    end

    if ifelse_with_time == true
        formula = _time_dependent_ifelse_to_bool(formula, model_SBML, ifelse_parameter_names)
    end
    return formula
end

function _get_side_activated_with_time(
        lhs_condition::String, rhs_condition::String, operator::String, time_in_rhs::Bool
    )::String
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

function _template_bool_picewise(
        bool_name::String, ifelse_arg1::String, ifelse_arg2::String, side_activated::String
    )::String
    activated = side_activated == "left" ? ifelse_arg1 : ifelse_arg2
    deactivated = side_activated == "left" ? ifelse_arg2 : ifelse_arg1
    formula = "((1 - 1" * bool_name * ") * (" * deactivated * ") + " *
        bool_name * "*(" * activated * "))"
    return formula
end

function _get_sign_time(formula::String)::Int64
    # Math expressions are stored in prefix notation (e.g. -(t, 5.0)). Meta.parse handles
    # both prefix and infix notation, which allows the direction to be inferred from the
    # expression tree
    expr = try
        Meta.parse(formula)
    catch
        nothing
    end
    if !(expr isa Expr && expr.head in [:error, :incomplete])
        sign_time = _get_sign_time(expr)
        !isnothing(sign_time) && return sign_time
    end

    # If a '-' does not appear in the formula the expression must increase with time. This
    # covers expressions the analysis above cannot handle (e.g. exp(t)). If a '-' appears
    # we might be able to infer direction, but infering in this situation is hard! - so
    # throw an error as the user should be able to write condition in a more easy manner
    # (avoid several sign changing minus signs)
    !occursin('-', formula) && return 1
    str_write = "For piecewise with time in condition we cannot infer direction for \
        $formula, that is if the condition value increases or decreases with time. This \
        happens if the formula contains a minus sign in the term where t appears."
    throw(SBMLSupport(str_write))
end

"""
    _get_sign_time(expr)::Union{Int64, Nothing}

Infer whether a Julia expression increases (`1`) or decreases (`-1`) with time.

If the direction cannot be inferred `nothing` is returned.
"""
function _get_sign_time(expr)::Union{Int64, Nothing}
    expr == :t && return 1
    !(expr isa Expr && expr.head == :call) && return nothing
    fn, args = expr.args[1], expr.args[2:end]
    if fn == :+
        return _get_sign_time_terms(args, fill(1, length(args)))
    elseif fn == :- && length(args) == 1
        sign_time = _get_sign_time(args[1])
        return isnothing(sign_time) ? nothing : -sign_time
    elseif fn == :- && length(args) == 2
        return _get_sign_time_terms(args, [1, -1])
    elseif fn == :*
        return _get_sign_time_factors(args)
    elseif fn == :/ && length(args) == 2
        # With time in the denominator the direction depends on the sign of the
        # denominator, which cannot be inferred
        _expr_has_time(args[2]) && return nothing
        return _get_sign_time_factors(args)
    end
    return nothing
end

function _get_sign_time_terms(args, signs)::Union{Int64, Nothing}
    sign_time = nothing
    for (arg, sign_term) in zip(args, signs)
        _expr_has_time(arg) == false && continue
        _sign_time = _get_sign_time(arg)
        isnothing(_sign_time) && return nothing
        _sign_time *= sign_term
        if isnothing(sign_time)
            sign_time = _sign_time
        elseif sign_time != _sign_time
            # Terms with time that change in opposite directions with time
            return nothing
        end
    end
    return sign_time
end

function _get_sign_time_factors(args)::Union{Int64, Nothing}
    itime = findall(_expr_has_time, args)
    length(itime) != 1 && return nothing
    sign_time = _get_sign_time(args[itime[1]])
    isnothing(sign_time) && return nothing
    # Direction can only be inferred if the sign of the remaining factors is known, which
    # is the case if they are numeric values
    for (i, arg) in pairs(args)
        i == itime[1] && continue
        value = _get_expr_value(arg)
        (isnothing(value) || value == 0) && return nothing
        sign_time *= value > 0 ? 1 : -1
    end
    return sign_time
end

function _get_expr_value(expr)::Union{Float64, Nothing}
    expr isa Number && return Float64(expr)
    if expr isa Expr && expr.head == :call && expr.args[1] == :- && length(expr.args) == 2
        value = _get_expr_value(expr.args[2])
        return isnothing(value) ? nothing : -value
    end
    return nothing
end

function _expr_has_time(expr)::Bool
    expr == :t && return true
    !(expr isa Expr) && return false
    return any(_expr_has_time, expr.args)
end

function _ifelse_to_event(id::String, condition::String, side_activated)::EventSBML
    assignments = [id * " = 1.0"]
    # When triggered we change the bool variable from 0 to 1. If side_activated = "right"
    # we activate the event when the ifelse condition goes from true to false. Reverse for
    # side_activated = "left". Thus, for side_activated = "right" we need to invert the
    # condition, as otherwise we mess up callback initialization where the bool variable
    # is set to 1 if the condition is true. It is also important to work with strict
    # inequality here to handle any edge-cases where the trigger time is t = 0
    if side_activated == "right"
        if any(occursin.(["<", "≤", "<="], condition))
            condition = replace(condition, r"≤|<=|<" => "≥")
        else
            condition = replace(condition, r"≥|>=|>" => "≤")
        end
    else
        if any(occursin.(["<", "≤", "<="], condition))
            condition = replace(condition, r"<=|<" => "≤")
        else
            condition = replace(condition, r">=|>" => "≥")
        end
    end
    event = EventSBML(
        id, condition, assignments, false, false, false, false, false, false, false, true
    )
    return event
end
