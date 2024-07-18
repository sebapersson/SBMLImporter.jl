"""
    parse_math(math_sbml::SBML.MathApply, libsbml_model::SBML.Model)::MathSBML

Parse SBML.jl math expression into a string with Julia syntax, and captures SBML IDs.

It is important for performance to capture the IDs of SBML species and parameters, as
these might be replaced later following the SBML specification. Therefore, the results
are stored in a `MathSBML` struct.

## Returns
- A `MathSBML` struct with the following relevant fields:
  - `formula::String`: The equations parsed from the SBML math.
  - `math_ids::Vector{String}`: Stores the math IDs.
  - `fns::Vector{String}`: Stores the functions applied to the expression
"""
function parse_math(math_sbml::SBML.MathApply, libsbml_model::SBML.Model)::MathSBML
    math_expression = MathSBML("", String[], String[], String[], false, false)
    formula = _parse_math(math_sbml, math_expression, libsbml_model)
    unique!(math_expression.math_idents)
    unique!(math_expression.fns)
    math_expression.formula = formula
    # For efficient insertion of functions into formulas the user defined functions
    # are stored. Note that inequality functions do not have a corresponding Julia version,
    # so they are treated as user functions
    for fn in math_expression.fns
        if haskey(libsbml_model.function_definitions, fn)
            push!(math_expression.user_fns, fn)
        end
    end
    for fn in ["gt", "lt", "geq", "leq"]
        push!(math_expression.user_fns, fn)
    end
    for ident in math_expression.math_idents
        if haskey(libsbml_model.reactions, ident)
            math_expression.has_reaction_ids = true
            break
        end
    end
    math_expression.has_rateOf = "rateOf" in math_expression.fns
    return math_expression
end
function parse_math(math_sbml::SBMLMathVariables, libsbml_model::SBML.Model)::MathSBML
    math_expression = MathSBML("", String[], String[], String[], false, false)
    formula = _parse_math(math_sbml, math_expression, libsbml_model)
    math_expression.formula = formula
    for ident in math_expression.math_idents
        if haskey(libsbml_model.reactions, ident)
            math_expression.has_reaction_ids = true
            break
        end
    end
    return math_expression
end
function parse_math(math_sbml::Nothing, libsbml_model::SBML.Model)::MathSBML
    return MathSBML("", String[], String[], String[], false, false)
end

function _parse_math(math_sbml::SBML.MathApply, math_expression::MathSBML,
                     libsbml_model::SBML.Model)::String
    fn = _parse_fn(math_sbml, libsbml_model.function_definitions)
    push!(math_expression.fns, fn)

    nargs = length(math_sbml.args)
    args = fill("", nargs)
    for i in 1:nargs
        args[i] = _parse_math(math_sbml.args[i], math_expression, libsbml_model)
    end

    # Varius edge cases
    # SBML MathML allows zero argument addition and multiplication.
    if fn in ["+", "*"] && nargs == 0
        return fn == "+" ? "0" : "1"
    end
    # Julia does not have a root function, hence for two args this is needed
    if fn == "sqrt" && nargs == 2
        return _apply(^, args[2], _apply(/, "1.0", args[1]))
    end
    # Gamma function requires a shift of +1; factorial(a) = Gamma(a+1)
    if fn == "SpecialFunctions.gamma"
        args[1] = _apply(+, args[1], "1")
    end
    # For gated logical funcitons 0 and 1 arguments are allowed. As these cases correspond
    # to different functions than the normal gated logicals, these are renamed
    if fn in ["and", "or", "if", "xor"] && nargs != 2
        fn *= string(nargs)
    end
    if fn == "piecewise" && nargs in [2, 4]
        @warn "Piecewise with $nargs (not 3) arguments is allowed but not recomended."
        fn *= string(nargs)
    end
    return fn * "(" * prod(args .* ", ")[1:(end - 2)] * ")"
end
_parse_math(math::SBML.MathVal, ::MathSBML, ::SBML.Model)::String = string(math.val)
_parse_math(math::SBML.MathTime, ::MathSBML, ::SBML.Model)::String = "t"
_parse_math(math::SBML.MathAvogadro, ::MathSBML, ::SBML.Model)::String = "6.02214179e23"
function _parse_math(math::SBML.MathIdent, math_expression::MathSBML, ::SBML.Model)::String
    id = string(math.id)
    push!(math_expression.math_idents, id)
    return id
end
function _parse_math(math::SBML.MathConst, ::MathSBML, ::SBML.Model)::String
    @assert math.id in ["pi", "exponentiale", "false", "true"] "Invalid arg constant $(arg.id) in SBML file"
    math.id == "exponentiale" && return "2.718281828459045"
    math.id == "pi" && return "3.1415926535897"
    return math.id
end

function _parse_fn(math_sbml::SBML.MathApply, sbml_functions::Dict)::String
    if haskey(sbml_functions, math_sbml.fn)
        return math_sbml.fn
    end

    @assert haskey(SBML_FN_INFO, math_sbml.fn) "$(math_sbml.fn) not among importable functions"
    fn, nallowed_args = SBML_FN_INFO[math_sbml.fn]
    nargs = length(math_sbml.args)
    if fn in ["delay", "rem", "implies", "div"]
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
