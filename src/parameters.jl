function parse_parameters!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (parameter_id, parameter) in libsbml_model.parameters
        model_SBML.parameters[parameter_id] = _parse_parameter(parameter_id, parameter)
    end
    return nothing
end

function _parse_parameter(parameter_id::String, parameter::SBML.Parameter)::ParameterSBML
    if parameter_id in FORBIDDEN_IDS
        throw(SBMLSupport("Parameter name $(parameter_id) is not allowed."))
    end
    formula = _parse_variable(parameter.value; default = "0.0")
    constant = _parse_bool(parameter.constant)
    initial_value = ""
    return ParameterSBML(parameter_id, constant, formula, initial_value,
        false, false, false, false, false, false)
end
