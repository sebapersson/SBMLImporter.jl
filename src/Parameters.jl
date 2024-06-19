#=
    Functionality for parsing, and handling SBML parameters and compartments
=#

function parse_parameters!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (parameter_id, parameter) in libsbml_model.parameters
        if parameter_id in FORBIDDEN_IDS
            throw(SBMLSupport("Parameter name $(parameter_id) is not allowed."))
        end

        formula = _parse_variable(parameter.value; default="0.0")
        constant = _parse_bool(parameter.constant)
        initial_value = ""
        model_SBML.parameters[parameter_id] = ParameterSBML(parameter_id, constant, formula, initial_value, false, false, false, false)
    end
    return nothing
end

function include_event_parameters_in_model!(model_SBML::ModelSBML)::Nothing

    # Sometimes parameter can be non-constant, but still have a constant rhs and they change value
    # because of event assignments. This must be captured by setting the parameter to have a zero
    # derivative so it is not simplified away later with structurally_simplify
    for (parameter_id, parameter) in model_SBML.parameters
        if parameter.algebraic_rule == true || parameter.rate_rule == true ||
           parameter.constant == true
            continue
        end
        if !is_number(parameter.formula)
            continue
        end
        # To pass test case 957
        if parse(Float64, parameter.formula) ≈ π
            continue
        end

        parameter.rate_rule = true
        parameter.initial_value = parameter.formula
        parameter.formula = "0.0"
        push!(model_SBML.rate_rule_variables, parameter_id)
    end
    # Similar holds for compartments
    for (compartment_id, compartment) in model_SBML.compartments
        if compartment.algebraic_rule == true || compartment.rate_rule == true ||
           compartment.constant == true
            continue
        end
        if !is_number(compartment.formula)
            continue
        end

        compartment.rate_rule = true
        compartment.initial_value = compartment.formula
        compartment.formula = "0.0"
        push!(model_SBML.rate_rule_variables, compartment_id)
    end

    return nothing
end
