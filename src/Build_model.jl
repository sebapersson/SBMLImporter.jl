"""
    build_SBML_model(libsbml_model::SBML.Model; ifelse_to_callback::Bool=true)::ModelSBML

Given the path to a SBML file, builds an intermediate SBML model struct.

The SBML model struct stores the information needed to create a ODESystem or ReactionSystem.

Rewriting ifelse to Boolean callbacks is strongly recommended (if possible).

For testing path_SBML can be the model as a string if model_as_string=true.
"""
function build_SBML_model(path_SBML::String; ifelse_to_callback::Bool = true,
                          model_as_string = true, inline_assignment_rules::Bool = true,
                          mass_action::Bool = false)::ModelSBML
    model_str = _get_model_as_str(path_SBML, model_as_string)
    check_support(path_SBML, model_as_string)

    # stoichiometryMath only occurs in SBML level 2. If stoichiometryMath occur in the
    # model it is converted to a level three model via SBML.jl libsbml functionality
    if occursin("stoichiometryMath", model_str) == false
        libsbml_model = readSBMLFromString(model_str)
    else
        libsbml_model = readSBMLFromString(model_str,
                                           doc -> begin
                                               set_level_and_version(3, 2)(doc)
                                               convert_promotelocals_expandfuns(doc)
                                           end)
    end

    return _build_SBML_model(libsbml_model, ifelse_to_callback, inline_assignment_rules,
                             mass_action)
end

function _build_SBML_model(libsbml_model::SBML.Model, ifelse_to_callback::Bool,
                           inline_assignment_rules::Bool, mass_action::Bool)::ModelSBML
    model_SBML = ModelSBML(libsbml_model)

    parse_species!(model_SBML, libsbml_model, mass_action)

    parse_parameters!(model_SBML, libsbml_model)

    parse_compartments!(model_SBML, libsbml_model)

    parse_functions!(model_SBML, libsbml_model)

    parse_rules!(model_SBML, libsbml_model)

    parse_events!(model_SBML, libsbml_model)

    parse_initial_assignments!(model_SBML, libsbml_model)

    parse_reactions!(model_SBML, libsbml_model)

    # Following the SBML standard reaction ids can appear in formulas, where they correspond
    # to the reaction kinetic_math expression
    replace_reactionid!(model_SBML)

    # Rewrite any time-dependent ifelse to boolean statements so that these can be expressed
    # as efficient DiscreteCallbacks later
    if ifelse_to_callback == true
        time_dependent_ifelse_to_bool!(model_SBML)
    end

    identify_algebraic_rule_variables!(model_SBML)

    add_conversion_factor_ode!(model_SBML, libsbml_model)

    # SBML allows inconstant compartment size, this must be adjusted if a specie is given in concentration
    # Must be after conversion factor, as the latter must be correctly handled in the transformation
    adjust_for_dynamic_compartment!(model_SBML)

    # Ensure that event participating parameters and compartments are not simplfied away when calling
    # structurally_simplify
    include_event_parameters_in_model!(model_SBML)

    # Per level3 rateOf can appear in any formula, and should be replaced with corresponding rate
    replace_rateOf!(model_SBML)

    # Inlining assignment rule variables makes the model less readable, however, if not done for
    # larger models Catalyst might crash due to a stack-overflow error
    if inline_assignment_rules == true
        inline_assignment_rules!(model_SBML)
    end

    # For ease of processing down the line storing all rule variables is convenient
    get_rule_variables!(model_SBML)

    return model_SBML
end
