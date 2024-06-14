function parse_functions!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (function_name, _function) in libsbml_model.function_definitions
        # Per SBML standard a function can be empty, in practice should be user misstake
        if isnothing(_function.body)
            @warn "Function $(function_name) does not have a function body and is ignored"
            continue
        end
        args = _get_function_args(_function)
        function_body = parse_math(_function.body.body, libsbml_model)
        model_SBML.functions[function_name] = FunctionSBML(args, function_body.formula)
    end

    # In SBML inequalities are provided as lt, gt, ..., they are basically SBML functions,
    # and are treated as such
    _add_inequality_functions!(model_SBML)
    return nothing
end

function _get_function_args(_function::SBML.FunctionDefinition)::Vector{String}
    isempty(_function.body.args) && return String[]
    return _function.body.args
end

function _add_inequality_functions!(model_SBML::ModelSBML)::Nothing
    # __x__ and __y__ to not end up using model ids
    model_SBML.functions["gt"] = FunctionSBML(["__x__", "__y__"], "__x__ > __y__")
    model_SBML.functions["lt"] = FunctionSBML(["__x__", "__y__"], "__x__ < __y__")
    model_SBML.functions["geq"] = FunctionSBML(["__x__", "__y__"], "__x__ >= __y__")
    model_SBML.functions["leq"] = FunctionSBML(["__x__", "__y__"], "__x__ <= __y__")
    return nothing
end

function insert_functions(formula::T, _functions::Dict)::T where {T <: AbstractString}
    # TODO: For early exit formula should track functions, can be done via math.jl,
    # and will make this function faster for larger models, as the haskey function
    # can be used here
    if !any(occursin.(keys(_functions) .* '(', formula))
        return formula
    end

    # TODO: As above, can be tracked with math.jl for faster early exit
    if replace_variable(formula, "and", "") != formula
        throw(SBMLSupport("and is not supported in SBML functions"))
    end
    if replace_variable(formula, "or", "") != formula
        throw(SBMLSupport("or is not supported in SBML functions"))
    end

    for (function_name, _function) in _functions
        # TODO: Early exit, to make more efficient as above let math.jl track
        # (will not need later as will loop over know functions in the expression)
        if !occursin(function_name * '(', formula)
            continue
        end

        nfunction_occurrences = _get_times_appear(function_name, formula)
        for _ in 1:nfunction_occurrences
            function_call = _extract_function_call(function_name, formula)
            function_insert = _get_expression_insert(function_call, _function)
            formula = replace(formula, function_call => function_insert; count=1)
        end
    end

    # Functions can be nested, handled via recursion
    if any(occursin.(keys(_functions), formula))
        formula = insert_functions(formula, _functions)
    end

    return formula
end

function _extract_function_call(name, formula)
    iname = findfirst(Regex("\\b" * name * "\\("), formula)
    iend = _find_indices_outside_paranthesis(')', formula[iname[1]:end])
    return formula[iname[1]:(iname[1]+iend[1]-1)]
end

function _get_expression_insert(function_call::String, _function::FunctionSBML)::String
    args_insert = _extract_args_insert(function_call)
    @unpack args, body = _function
    @assert length(args) == length(args_insert) "Number of arguments to insert does not match SBML function"
    function_inserted = deepcopy(body)
    for i in eachindex(args_insert)
        function_inserted = replace_variable(function_inserted, args[i], args_insert[i])
    end
    return function_inserted
end

function _extract_args_insert(function_call::String)::Vector{String}
    # Args are , separated, on form f(arg1, arg2, ..., argn), which is used for extracting
    icomma = _find_indices_outside_paranthesis(',', function_call; start_depth=-1)
    istart_args = findfirst('(', function_call) + 1
    args_insert = _split_by_indicies(function_call, icomma; istart=istart_args, iend=1)
    # To ensure correct function evaluation order as args can be expression, paranthesis is
    # added around each arg. E.g. test case 1770
    args_insert = '(' .* replace.(args_insert, " " => "") .* ')'
    return args_insert
end
