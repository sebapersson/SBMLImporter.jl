function parse_compartments!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (compartment_id, compartment) in libsbml_model.compartments
        if compartment_id in FORBIDDEN_IDS
            throw(SBMLSupport("Compartment name $(parameter_id) is not allowed."))
        end

        size = _parse_variable(compartment.size; default="1.0")
        constant = _parse_bool(compartment.constant)
        initial_value = ""
        model_SBML.compartments[compartment_id] = CompartmentSBML(compartment_id, constant, size, initial_value, false, false, false, false, false)
    end
    return nothing
end
