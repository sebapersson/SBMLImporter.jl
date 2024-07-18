function _apply(operator::Function, x::String, y::String)::String
    @assert operator in [+, -, *, /, ^] "$operator is not allowed for building math expression"
    return string(operator) * '(' * x * ", " * y * ')'
end

function _template_value_map(id::String, value::String)::String
    return "\t" * id * " =>" * value * ",\n"
end

function _template_rate_rule(id::String, formula::String)::String
    return "\t\tD(" * id * ") ~ " * formula * ",\n"
end

function _template_assignment_rule(id::String, formula::String)::String
    return "\t\t" * id * " ~ " * formula * ",\n"
end

function _template_reaction(reactants::String, products::String, r_S::String, p_S::String,
                            propensity::String, is_massaction::Bool)::String
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
    return _apply(*, s, c_scaling)
end

function _template_ode_reaction(s::String, c_scaling::String, propensity::String,
                                which_side::Symbol)::String
    s == "nothing" && return ""
    @assert which_side in [:reactant, :product] "$(which_side) is an invalid reaction side"
    sign = which_side == :product ? '+' : '-'
    # c_scaling either "" or "/...", hence the 1
    p1 = _apply(*, c_scaling, propensity)
    p2 = _apply(*, "1", p1)
    return _apply(*, sign * s, p2)
end

function _template_tstops(tstops::Vector{String})::String
    if isempty(tstops)
        tstops_vec = "Float64[]"
    else
        tstops_vec = "Float64[" .* prod([t * ", " for t in tstops]) * "]"
    end
    return tstops_vec
end

function _template_condition(condition::String, discrete_callback::Bool,
                             name::String)::String
    # Discrete callbacks are only activated when condition transition from false to true
    # Continious callbacks are activated via root-finding (hence the -)
    if discrete_callback == true
        condition_f = "\nfunction _condition_" * name * "(u, t, integrator, from_neg)\n"
        condition_bool = condition
        condition = "\tcond = " * condition_bool * " && from_neg[1] == true\n"
        condition *= "\n\t\tfrom_neg[1] = !(" * condition_bool * ")\n"
        condition *= "\t\treturn cond"
    elseif discrete_callback == false
        condition_f = "\nfunction condition_" * name * "(u, t, integrator)\n"
        condition = replace(condition, r"<=|>=|>|<|≤|≥" => "-")
    end
    condition_f *= condition * "\nend"
    return condition_f
end

function _template_affect(event::EventSBML, specie_ids::Vector{String},
                          parameter_ids::Vector{String}; only_body::Bool = false)::String
    if only_body == false
        affect_f = "function affect_" * event.name * "!(integrator)\n"
    else
        affect_f = ""
    end
    affect_f *= "\tutmp = similar(integrator.u)\n\tutmp .= integrator.u\n"

    for affect in event.formulas
        affect_lhs, affect_rhs = string.(split(affect, "="))
        affect_lhs = _ids_to_callback_syntax(affect_lhs, specie_ids, :specie)
        affect_rhs = _ids_to_callback_syntax(affect_rhs, specie_ids, :specie;
                                             integrator = false, utmp = true)
        affect_eq = affect_lhs * "=" * affect_rhs
        affect_eq = _ids_to_callback_syntax(affect_eq, parameter_ids, :parameter)
        affect_f *= "\t" * affect_eq * "\n"
    end

    # Needed for events to work with JumpProblem for Gillespie type simulations
    if only_body == false
        affect_f *= "\tif integrator.sol.prob isa SciMLBase.DiscreteProblem\n"
        affect_f *= "\t\treset_aggregated_jumps!(integrator)\n\tend\n"
        affect_f *= "end"
    end
    return affect_f
end

function _template_init(event::EventSBML, condition::String, affect_body::String,
                        tstops::String, first_callback::Bool,
                        discrete_callback::Bool)::String
    init = "function init_" * event.name * "(c,u,t,integrator)\n"
    # For DiscreteCallback's we might need to support tstops, these can be computed in the
    # init
    if first_callback == true
        init *= "\t_tstops = " * tstops * "\n"
        init *= "tstops = _tstops[@.((integrator.tdir * _tstops > integrator.tdir * integrator.sol.prob.tspan[1])*(integrator.tdir *_tstops < integrator.tdir * integrator.sol.prob.tspan[2]))]\n"
        init *= "\ttstops = isempty(tstops) ? tstops : vcat(minimum(tstops) / 2.0, tstops)\n"
        init *= "\tSciMLBase.add_tstop!.((integrator,), tstops)"
    end

    # Note, init for events can only be triggered if the SBML event has
    # trigger_initial_value (which is encoded by default for SBML piecewise)
    skip_init = (discrete_callback && event.is_ifelse == false)
    if event.trigger_initial_value == true || skip_init
        init *= "\nend"
        return init
    end

    init *= "\n\tcond = " * condition * "\n"
    init *= "\tif cond == true\n"
    init *= "\t" * affect_body * "\n\tend"
    init *= "\nend"
    return init
end

function _template_conc_dynamic_c(dndt::String, V::String, n_id::String,
                                  dVdt::String)::String
    # Math expressions are built via, e.g., +(1, 2) to avoid paranthesis problems
    V2 = _apply(^, V, "2")
    p1 = _apply(/, dndt, V)
    p2 = _apply(*, _apply(/, n_id, V2), dVdt)
    return _apply(-, p1, p2)
end

function _template_amount_dynamic_c(dcdt::String, V::String, specie_id::String,
                                    dVdt::String)::String
    p1 = _apply(*, dcdt, V)
    p2 = _apply(*, specie_id, _apply(/, dVdt, V))
    return _apply(+, p1, p2)
end
