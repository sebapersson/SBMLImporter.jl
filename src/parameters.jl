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
