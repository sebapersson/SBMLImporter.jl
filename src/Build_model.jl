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
    if model_as_string == false
        if !isfile(path_SBML)
            throw(SBMLSupport("$path_SBML is not the path to a SBML file"))
        end
        f = open(path_SBML, "r")
            model_str = read(f, String)
        close(f)
    else
        model_str = path_SBML
    end

    if has_event_with_delay(model_str) == true
        throw(SBMLSupport("Events with delay are not supported"))
    end
    if has_event_with_priority(model_str) == true
        throw(SBMLSupport("Events with priority are not supported"))
    end
    if has_fast_reaction(model_str) == true
        throw(SBMLSupport("Fast reactions are not supported"))
    end
    if is_hierarchical(model_str) == true
        throw(SBMLSupport("Hierarchical models are not supported"))
    end
    if is_fba(model_str) == true
        throw(SBMLSupport("FBA models are not supported. Checkout COBREXA.jl"))
    end

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
    conversion_factor = isnothing(libsbml_model.conversion_factor) ? "" :
                        libsbml_model.conversion_factor

    # Convert model name to a valid Julia string
    model_name = isnothing(libsbml_model.name) ? "SBML_model" : libsbml_model.name
    model_name = replace(model_name, " " => "_")

    # Specie reference ids can sometimes appear in math expressions, then they should be replaced
    # by the stoichometry for corresponding reference id, searching for specie reference ids
    # can be computationally demanding if the list of such ids is rebuilt, thus here the
    # list is built only once in the beginning
    specie_reference_ids = get_specie_reference_ids(libsbml_model)

    # An intermedidate struct storing relevant model informaiton needed for
    # formulating a ReactionSystem and callback functions
    model_SBML = ModelSBML(model_name; specie_reference_ids = specie_reference_ids,
                           conversion_factor = conversion_factor)

    parse_species!(model_SBML, libsbml_model, mass_action)

    parse_parameters!(model_SBML, libsbml_model)

    parse_compartments!(model_SBML, libsbml_model)

    parse_functions!(model_SBML, libsbml_model)

    parse_rules!(model_SBML, libsbml_model)

    parse_SBML_events!(model_SBML, libsbml_model)

    # Positioned after rules since some assignments may include functions
    parse_SBML_initial_assignments!(model_SBML, libsbml_model)

    parse_SBML_reactions!(model_SBML, libsbml_model)

    # Given the SBML standard reaction id can sometimes appear in the reaction
    # formulas, here the correpsonding id is overwritten with a math expression
    replace_reactionid!(model_SBML)

    # Rewrite any time-dependent ifelse to boolean statements such that we can express these as events.
    # This is recomended, as it often increases the stabillity when solving the ODE, and decreases run-time
    if ifelse_to_callback == true
        time_dependent_ifelse_to_bool!(model_SBML)
    end

    identify_algebraic_rule_variables!(model_SBML)

    adjust_conversion_factor!(model_SBML, libsbml_model)

    # SBML allows inconstant compartment size, this must be adjusted if a specie is given in concentration
    # Must be after conversion factor, as the latter must be correctly handled in the transformation
    adjust_for_dynamic_compartment!(model_SBML)

    # Ensure that event participating parameters and compartments are not simplfied away when calling
    # structurally_simplify
    include_event_parameters_in_model!(model_SBML)

    # Per level3 rateOf can appear in any formula, and should be replaced with corresponding rate
    replace_rateOf!(model_SBML)

    # Up to this point for models parameter the SBML rem or div function might appear in the model dynamics,
    # these are non differentialble, discrete and not allowed
    has_rem_or_div(model_SBML)

    # Inlining assignment rule variables makes the model less readable, however, if not done for
    # larger models Catalyst might crash due to a stack-overflow error
    if inline_assignment_rules == true
        inline_assignment_rules!(model_SBML)
    end

    # For ease of processing down the line storing all rule variables is convenient
    get_rule_variables!(model_SBML)

    return model_SBML
end
