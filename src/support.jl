function check_support(model_str::String)::Nothing
    if _has_event_with_delay(model_str) == true
        throw(SBMLSupport("Events with delay are not supported"))
    end
    if _has_event_with_priority(model_str) == true
        throw(SBMLSupport("Events with priority are not supported"))
    end
    if _has_fast_reaction(model_str) == true
        throw(SBMLSupport("Fast reactions are not supported"))
    end
    if _is_hierarchical(model_str) == true
        throw(SBMLSupport("Hierarchical models are not supported"))
    end
    if _is_fba(model_str) == true
        throw(SBMLSupport("FBA models are not supported. Checkout COBREXA.jl"))
    end

    return nothing
end

function _has_event_with_delay(file_text::String)::Bool
    istart = findall("<listOfEvents>", file_text)
    iend = findall("</listOfEvents>", file_text)
    if isnothing(istart) || isempty(istart)
        return false
    end

    return occursin("delay", file_text[istart[1][1]:iend[1][end]])
end

function _has_event_with_priority(file_text::String)::Bool
    istart = findall("<listOfEvents>", file_text)
    iend = findall("</listOfEvents>", file_text)
    if isnothing(istart) || isempty(istart)
        return false
    end
    return occursin("priority", file_text[istart[1][1]:iend[1][end]])
end

function _has_fast_reaction(file_text::String)::Bool
    return occursin(r"fast=\"true\"", file_text)
end

function _is_hierarchical(file_text::String)::Bool
    return occursin("comp:", file_text)
end

function _is_fba(file_text::String)::Bool
    return occursin("fbc:", file_text)
end
