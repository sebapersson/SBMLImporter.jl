"""
    build_SBML_model(libsbml_model::SBML.Model; ifelse_to_callback::Bool=true)::ModelSBML

Given the path to a SBML file build an intermediate SBML model struct.

The SBML model struct stores the information needed to create a ODESystem or ReactionSystem

Rewriting ifelse to Boolean callbacks is strongly recomended if possible.

For testing path_SBML can be the model as a string if model_as_string=true
"""
function build_SBML_model(path_SBML::String; ifelse_to_callback::Bool=true, model_as_string=true)::ModelSBML

    if model_as_string == false
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
        libsbml_model = readSBMLFromString(model_str, doc -> begin
                                 set_level_and_version(3, 2)(doc)
                                 convert_promotelocals_expandfuns(doc)
                                 end)
    end

    return _build_SBML_model(libsbml_model, ifelse_to_callback)
end



function _build_SBML_model(libsbml_model::SBML.Model, ifelse_to_callback::Bool)::ModelSBML

    # An intermedidate struct storing relevant model informaiton needed for
    # formulating an ODESystem and callback functions
    conversion_factor = isnothing(libsbml_model.conversion_factor) ? "" : libsbml_model.conversion_factor
    model_name = isnothing(libsbml_model.name) ? "SBML_model" : libsbml_model.name
    model_name = replace(model_name, " " => "_")
    model_SBML = ModelSBML(model_name,
                           Dict{String, SpecieSBML}(),
                           Dict{String, ParameterSBML}(),
                           Dict{String, CompartmentSBML}(),
                           Dict{String, EventSBML}(),
                           Dict{String, ReactionSBML}(),
                           Dict{String, Vector{String}}(), # SBML reactions
                           Dict{String, String}(), # Algebraic rules
                           Dict{String, String}(), # Generated id:s
                           Dict{String, String}(), # Piecewise to ifelse_expressions
                           Dict{String, String}(), # Ifelse to bool expression
                           Dict{String, Vector{String}}(), # Ifelse parameters
                           Vector{String}(undef, 0), # Rate rule variables
                           Vector{String}(undef, 0), # Assignment rule variables
                           Vector{String}(undef, 0), # Algebraic rule variables
                           Vector{String}(undef, 0), # Species_appearing in reactions
                           Vector{String}(undef, 0), 
                           conversion_factor) # Variables with piecewise

    parse_SBML_species!(model_SBML, libsbml_model)

    parse_SBML_parameters!(model_SBML, libsbml_model)

    parse_SBML_compartments!(model_SBML, libsbml_model)

    parse_SBML_functions!(model_SBML, libsbml_model)

    parse_SBML_rules!(model_SBML, libsbml_model)

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

    return model_SBML
end