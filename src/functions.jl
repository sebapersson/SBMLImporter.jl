function parse_functions!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (function_name, _function) in libsbml_model.function_definitions
        # Per SBML standard a function can be empty, in practice should be user misstake
        if isnothing(_function.body)
            @warn "Function $(function_name) does not have a function body and is ignored"
            continue
        end
        args = _get_function_args(_function)
        function_body = parse_math(_function.body.body, libsbml_model)
        # To avoid problems with inserting a function into a math expression the function
        # args are not allowed to match any other SBML model variables. Hence, args are
        # renamed to have underscore names
        for (i, arg) in pairs(args)
            args[i] = "__" * arg * "__"
            function_body.formula = _replace_variable(function_body.formula, arg, args[i])
        end
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
    model_SBML.functions["geq"] = FunctionSBML(["__x__", "__y__"], "__x__ ≥ __y__")
    model_SBML.functions["leq"] = FunctionSBML(["__x__", "__y__"], "__x__ ≤ __y__")
    return nothing
end

function insert_functions(formula::T, _functions::Dict,
        fns_reaplce::Vector{String})::T where {T <: AbstractString}
    isempty(fns_reaplce) && return formula

    nfunctions_replaced::Int64 = 0
    for function_name in fns_reaplce
        # As (due to nesting) functions are pased via recursion a function_name might not
        # appear in the formula
        !occursin(function_name * '(', formula) && continue

        _function = _functions[function_name]
        nfunction_occurrences = _get_times_appear(function_name, formula)
        for _ in 1:nfunction_occurrences
            function_call = _extract_function_call(function_name, formula)
            function_insert = _get_expression_insert(function_call, _function)
            formula = replace(formula, function_call => function_insert; count = 1)
            nfunctions_replaced += 1
        end
    end
    # For fast exit when hitting recrusion bottom
    nfunctions_replaced == 0 && return formula
    # Functions can be nested, with a function returning a function. In this case,
    # recursively all possible functions in the model
    if any(occursin.(keys(_functions), formula))
        formula = insert_functions(formula, _functions, collect(keys(_functions)))
    end
    return formula
end

function _extract_function_calls(name, formula)
    inames = findall(Regex("\\b" * name * "\\("), formula)
    out = fill("", length(inames))
    for i in eachindex(out)
        iname = inames[i]
        iend = _find_indices_outside_paranthesis(')', formula[iname[1]:end])
        out[i] = formula[iname[1]:(iname[1] + iend[1] - 1)]
    end
    return out
end

function _extract_function_call(name, formula)
    iname = findfirst(Regex("\\b" * name * "\\("), formula)
    iend = _find_indices_outside_paranthesis(')', formula[iname[1]:end])
    return formula[iname[1]:(iname[1] + iend[1] - 1)]
end

function _get_expression_insert(function_call::String, _function::FunctionSBML)::String
    args_insert = _extract_args_insert(function_call)
    @unpack args, body = _function
    @assert length(args)==length(args_insert) "Number of arguments to insert does not match SBML function"
    function_inserted = deepcopy(body)
    for i in eachindex(args_insert)
        function_inserted = _replace_variable(function_inserted, args[i], args_insert[i])
    end
    return function_inserted
end

function _extract_args_insert(function_call::String)::Vector{String}
    # Args are , separated, on form f(arg1, arg2, ..., argn), which is used for extracting
    icomma = _find_indices_outside_paranthesis(',', function_call; start_depth = -1)
    istart_args = findfirst('(', function_call) + 1
    args_insert = _split_by_indicies(function_call, icomma; istart = istart_args, iend = 1)
    return args_insert
end
