function _template_value_map(id::String, value::String)::String
    return "\t" * id * " =>" * value * ",\n"
end

function _template_rate_rule(id::String, formula::String)::String
    return "\t\tD(" * id * ") ~ " * formula * ",\n"
end

function _template_assignment_rule(id::String, formula::String)::String
    return "\t\t" * id * " ~ " * formula * ",\n"
end

function _template_reaction(reactants::String, products::String, r_S::String, p_S::String, propensity::String, is_massaction::Bool)::String
    reaction = "\t\t"
    if is_massaction
        reaction *= "SBMLImporter.update_rate_reaction("
    end
    reaction *= "Catalyst.Reaction(" *
        propensity * ", " *
        reactants * ", " *
        products * ", " *
        r_S * ", " *
        p_S
    if is_massaction
        reaction *= "; only_use_rate=false)),\n"
    else
        reaction *= "; only_use_rate=true),\n"
    end
    return reaction
end

function _template_stoichiometry(s::String, c_scaling::String)::String
    s == "nothing" && return "nothing"
    return s * c_scaling
end

function _template_ode_reaction(s::String, c_scaling::String, propensity::String, which_side::Symbol)::String
    s == "nothing" && return ""
    @assert which_side in [:reactant, :product] "$(which_side) is an invalid reaction side"
    sign = which_side == :product ? '+' : '-'
    # c_scaling either "" or "/...", hence the 1
    return sign * s * "*1" * c_scaling * "*(" * propensity * ")"
end
