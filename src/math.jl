"""
    parse_math(math_sbml::SBML.MathApply, libsbml_model::SBML.Model)::MathSBML

Parse SBML.jl math expression into a string with Julia syntax, and captures SBML IDs.

It is important for performance to capture the IDs of SBML species and parameters, as
these might be replaced later following the SBML specification. Therefore, the results
are stored in a `MathSBML` struct that during parsing (which may involve recursion),
this struct tracks math IDs.

## Returns
- A `MathSBML` struct with the following relevant fields:
  - `formula::String`: The equations parsed from the SBML math.
  - `math_ids::Vector{String}`: Stores the math IDs.
"""
function parse_math(math_sbml::SBML.MathApply, libsbml_model::SBML.Model)::MathSBML
    math_expression = MathSBML("", String[], "", String[])
    _parse_math!(math_expression, math_sbml, libsbml_model)
    unique!(math_expression.math_idents)
    return math_expression
end
function parse_math(math_sbml::SBMLMathVariables, libsbml_model::SBML.Model)::MathSBML
    formula = _parse_arg(math_sbml)
    return MathSBML(formula, String[formula], "", String[])
end
function parse_math(math_sbml::Nothing, libsbml_model::SBML.Model)::MathSBML
    return MathSBML("", String[], "", String[])
end

function _parse_math!(math_expression::MathSBML, math_sbml::SBML.MathApply, libsbml_model::SBML.Model; parse_fn_arg::Bool=false)::Nothing
    fn = _parse_fn(math_sbml, libsbml_model.function_definitions)
    _parse_args!(math_expression, math_sbml.args, libsbml_model)
    args = math_expression.args
    nargs = length(args)

    # Varius edge cases
    # SBML MathML allows zero argument addition and multiplication.
    if fn in ["+", "*"] && nargs == 0
        _formula = fn == "+" ? "0" : "1"
        _update_formula!(math_expression, _formula, parse_fn_arg)
        return nothing
    end
    # Must Julia functions are of type fn(args...) except +, *, -, /, ^.
    if fn in ["+", "-", "*", "/", "^"] && nargs == 2
        _formula = "(" * args[1] * ")" * fn * "(" * args[2] * ")"
        _update_formula!(math_expression, _formula, parse_fn_arg)
        return nothing
    end
    # Julia does not have a root function, hence for two args this is needed
    if fn == "sqrt" && nargs == 2
        _formula = "(" * args[2] * ")^(1.0 / " * args[1] * ")"
        _update_formula!(math_expression, _formula, parse_fn_arg)
        return nothing
    end
    # Gamma function requires a shift of +1; factorial(a) = Gamma(a+1)
    if fn == "SpecialFunctions.gamma"
        args[1] *= " + 1"
    end

    _formula = fn * "(" * prod(args .* ", ")[1:end-2] * ")"
    _update_formula!(math_expression, _formula, parse_fn_arg)
    return nothing
end

function _parse_fn(math_sbml::SBML.MathApply, sbml_functions::Dict)::String
    if haskey(sbml_functions, math_sbml.fn)
        return math_sbml.fn
    end

    @assert haskey(SBML_FN_INFO, math_sbml.fn) "$(math_sbml.fn) not among importable functions"
    fn, nallowed_args = SBML_FN_INFO[math_sbml.fn]
    nargs = length(math_sbml.args)
    if fn in ["delay", "rem", "implies"]
        throw(SBMLSupport("SBML models with function $fn are not supported"))
    end
    if fn in ["and", "xor", "or", "lt", "leq", "gt", "geq", "eq"] && nargs > 2
        throw(SBMLSupport("$fn with three conditions is not supported"))
    end
    if fn == "SpecialFunctions.gamma"
        @warn "Factorial in the ODE model. SBMLImporter can handle factorials, but, " *
              "solving the ODEs with factorial is numerically hard."
    end
    # Only empty for functions with arbibrary args
    if !isempty(nallowed_args)
        @assert nargs in nallowed_args "Function $fn with $nargs incorrect numbers of args"
    end
    return fn
end

function _parse_args!(math_expression::MathSBML, args_sbml, libsbml_model::SBML.Model)::Nothing
    args = Vector{String}(undef, length(args_sbml))
    for (i, _arg) in pairs(args_sbml)
        # SBML.MathApply triggers a recursion and is treated differently
        if _arg isa SBML.MathApply
            _parse_arg!(math_expression, _arg, libsbml_model)
            args[i] = math_expression.tmp_arg
        else
            args[i] = _parse_arg(_arg)
        end
        if _arg isa SBML.MathIdent
            push!(math_expression.math_idents, args[i])
        end
    end
    math_expression.args = args
    return nothing
end

_parse_arg(arg::SBML.MathVal)::String = string(arg.val)
_parse_arg(arg::SBML.MathIdent)::String = string(arg.id)
_parse_arg(arg::SBML.MathTime)::String = "t"
_parse_arg(arg::SBML.MathAvogadro)::String = "6.02214179e23"
function _parse_arg(arg::SBML.MathConst)::String
   @assert arg.id in ["pi", "exponentiale", "false", "true"] "Invalid arg constant $(arg.id) in SBML file"
    arg.id == "exponentiale" && return "2.718281828459045"
    arg.id == "pi" && return "3.1415926535897"
    return arg.id
end

function _parse_arg!(math_expression::MathSBML, arg::SBML.MathApply, libsbml_model::SBML.Model)::Nothing
    _parse_math!(math_expression, arg, libsbml_model; parse_fn_arg=true)
end

function _update_formula!(math_expression::MathSBML, formula::String, parse_fn_arg::Bool)
    if parse_fn_arg == false
        math_expression.formula = formula
    end
    if parse_fn_arg == true
        math_expression.tmp_arg = formula
    end
end
