"""
    parse_SBML(path, massaction::Bool; kwargs...)::ModelSBML

Parse the model into an intermediate SBML model struct that stores all information needed
to create a ReactionSystem.

For testing path can be the model as a string if model_as_string=true.
"""
function parse_SBML(
        path::String, massaction::Bool; ifelse_to_callback::Bool = true,
        model_as_string = true, inline_assignment_rules::Bool = true,
        inline_kineticlaw_parameters::Bool = true
    )::ModelSBML
    model_str = _get_model_as_str(path, model_as_string)
    check_support(path)

    # stoichiometryMath only occurs in SBML level 2. If stoichiometryMath occur in the
    # model it is converted to a level three model via SBML.jl libsbml functionality
    if occursin("stoichiometryMath", model_str) == false
        libsbml_model = readSBMLFromString(model_str)
    else
        libsbml_model = readSBMLFromString(
            model_str,
            doc -> begin
                SBML.set_level_and_version(3, 2)(doc)
                SBML.convert_promotelocals_expandfuns(doc)
            end
        )
    end

    return _parse_SBML(
        libsbml_model, ifelse_to_callback, inline_assignment_rules,
        inline_kineticlaw_parameters, massaction
    )
end

function _parse_SBML(
        libsbml_model::SBML.Model, ifelse_to_callback::Bool, inline_assignment_rules::Bool,
        inline_kineticlaw_parameters::Bool, massaction::Bool
    )::ModelSBML
    model_SBML = ModelSBML(libsbml_model)

    parse_species!(model_SBML, libsbml_model, massaction)

    parse_parameters!(model_SBML, libsbml_model)

    parse_compartments!(model_SBML, libsbml_model)

    parse_functions!(model_SBML, libsbml_model)

    parse_rules!(model_SBML, libsbml_model)

    parse_events!(model_SBML, libsbml_model)

    parse_initial_assignments!(model_SBML, libsbml_model)

    parse_reactions!(model_SBML, libsbml_model, inline_kineticlaw_parameters)

    # Following the SBML standard reaction ids can appear in formulas, where they correspond
    # to the reaction kinetic_math expression
    replace_ident!(model_SBML, libsbml_model, :reactionid)

    # Rewrite any time-dependent ifelse to boolean statements so that these can be expressed
    # as efficient DiscreteCallbacks later
    if ifelse_to_callback == true
        time_dependent_ifelse_to_bool!(model_SBML)
    end

    identify_algebraic_rule_variables!(model_SBML)

    add_conversion_factor_ode!(model_SBML, libsbml_model)

    # Ensure that event participating parameters and compartments are not simplified away
    # with structurally_simplify
    force_include_event_variables!(model_SBML)

    # SBML allows inconstant compartment size, this must be adjusted if a specie is given
    # in concentration. Must be after conversion factor, as the latter must be correctly
    # handled in the transformation
    adjust_for_dynamic_compartment!(model_SBML)

    # Per level3 rateOf can appear in any formula, and should be replaced with corresponding
    # rate
    replace_ident!(model_SBML, libsbml_model, :rateOf)

    # Specie reference ids can appear in formula, and need to be replaced
    replace_ident!(model_SBML, libsbml_model, :specieref)

    # Inlining assignment rule variables makes the model less readable, however, if not
    # done for larger models Catalyst might crash due to a stack-overflow error
    if inline_assignment_rules == true
        inline_assignment_rules!(model_SBML)
    end

    # For ease of processing down the line storing all rule variables is convenient
    get_rule_variables!(model_SBML)

    return model_SBML
end
