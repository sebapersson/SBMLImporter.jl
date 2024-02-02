function has_event_with_delay(file_text::String)::Bool
    istart = findall("<listOfEvents>", file_text)
    iend = findall("</listOfEvents>", file_text)
    if isnothing(istart) || isempty(istart)
        return false
    end

    @assert length(istart) == 1 && length(iend) == 1

    if occursin("delay", file_text[istart[1][1]:iend[1][end]])
        return true
    else
        return false
    end
end

function has_event_with_priority(file_text::String)::Bool
    istart = findall("<listOfEvents>", file_text)
    iend = findall("</listOfEvents>", file_text)
    if isnothing(istart) || isempty(istart)
        return false
    end

    @assert length(istart) == 1 && length(iend) == 1

    if occursin("priority", file_text[istart[1][1]:iend[1][end]])
        return true
    else
        return false
    end
end

function is_hierarchical(file_text::String)::Bool
    if occursin("comp:", file_text)
        return true
    else
        return false
    end
end

function is_fba(file_text::String)::Bool
    if occursin("fbc:", file_text)
        return true
    else
        return false
    end
end

function has_fast_reaction(file_text::String)::Bool
    if occursin(r"fast=\"true\"", file_text)
        return true
    end
    return false
end

function has_rem_or_div(model_SBML::ModelSBML)::Bool
    for (parameter_id, parameter) in model_SBML.parameters
        if occursin("div", parameter.formula)
            throw(SBMLSupport("quotient function in SBML model is not supported"))
        end
        if occursin("div", parameter.formula)
            throw(SBMLSupport("quotient function in SBML model is not supported"))
        end
    end
    return false
end
