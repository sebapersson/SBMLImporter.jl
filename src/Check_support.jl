function has_event_with_delay(path_SBML::String)::Bool

    f = open(path_SBML, "r")
        file_text = read(f, String)
    close(f)

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


function has_event_with_priority(path_SBML::String)::Bool

    f = open(path_SBML, "r")
        file_text = read(f, String)
    close(f)

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


function is_hierarchical(path_SBML::String)::Bool

    f = open(path_SBML, "r")
        file_text = read(f, String)
    close(f)

    if occursin("comp:", file_text)
        return true
    else
        return false
    end
end


function is_fba(path_SBML::String)::Bool

    f = open(path_SBML, "r")
        file_text = read(f, String)
    close(f)

    if occursin("fbc:", file_text)
        return true
    else
        return false
    end
end



function has_fast_reaction(path_SBML::String)::Bool

    f = open(path_SBML, "r")
        file_text = read(f, String)
    close(f)

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