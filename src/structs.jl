mutable struct SpecieSBML
    const name::String
    const boundary_condition::Bool
    const constant::Bool
    initial_value::String
    formula::String
    const compartment::String
    const conversion_factor::String
    const unit::Symbol
    const only_substance_units::Bool
    assignment_rule::Bool
    rate_rule::Bool
    algebraic_rule::Bool
    has_reaction_ids::Bool
    has_rateOf::Bool
    has_specieref::Bool
end

mutable struct ParameterSBML
    const name::String
    const constant::Bool
    formula::String
    initial_value::String
    assignment_rule::Bool
    rate_rule::Bool
    algebraic_rule::Bool
    has_reaction_ids::Bool
    has_rateOf::Bool
    has_specieref::Bool
end

mutable struct CompartmentSBML
    const name::String
    constant::Bool
    formula::String
    initial_value::String
    assignment_rule::Bool
    rate_rule::Bool
    algebraic_rule::Bool
    has_reaction_ids::Bool
    has_rateOf::Bool
    has_specieref::Bool
end

struct FunctionSBML
    args::Vector{String}
    body::String
end

mutable struct EventSBML
    const name::String
    trigger::String
    const formulas::Vector{String}
    const trigger_initial_value::Bool
    const has_reaction_ids_trigger::Bool
    const has_reaction_ids_assignments::Bool
    const has_rateOf_trigger::Bool
    const has_rateOf_assignments::Bool
    const has_specieref_trigger::Bool
    const has_specieref_assignments::Bool
    const is_ifelse::Bool
end

mutable struct ReactionSBML
    const name::String
    kinetic_math::String
    const products::Vector{String}
    const products_stoichiometry::Vector{String}
    const reactants::Vector{String}
    const reactants_stoichiometry::Vector{String}
    const stoichiometry_mass_action::Bool
    const has_assignment_rule_variable::Bool
    const has_reaction_ids::Bool
    const has_rateOf::Bool
    const has_specieref::Bool
end

mutable struct AlgebraicRuleSBML
    formula::String
    const idents::Vector{String}
    const has_rateOf::Bool
    const has_specieref::Bool
end

struct ModelSBML
    name::String
    species::Dict{String, SpecieSBML}
    parameters::Dict{String, ParameterSBML}
    compartments::Dict{String, CompartmentSBML}
    events::Dict{String, EventSBML}
    reactions::Dict{String, ReactionSBML}
    functions::Dict{String, FunctionSBML}
    algebraic_rules::Dict{String, AlgebraicRuleSBML}
    generated_ids::Dict{String, String}
    piecewise_expressions::Dict{String, String}
    ifelse_bool_expressions::Dict{String, String}
    rate_rule_variables::Vector{String}
    assignment_rule_variables::Vector{String}
    algebraic_rule_variables::Vector{String}
    species_in_reactions::Vector{String}
    variables_with_piecewise::Vector{String}
    conversion_factor::String
    specie_reference_ids::Vector{String}
    rule_variables::Vector{String}
    libsbml_rule_variables::Vector{String}
    reactionids::Vector{String}
end
function ModelSBML(name::String; specie_reference_ids::Vector{String} = String[],
                   conversion_factor::String = "", libsbml_rule_variables::Vector{String} = String[],
                   reactionids::Vector{String} = String[])::ModelSBML
    model_SBML = ModelSBML(name,
                           Dict{String, SpecieSBML}(),
                           Dict{String, ParameterSBML}(),
                           Dict{String, CompartmentSBML}(),
                           Dict{String, EventSBML}(),
                           Dict{String, ReactionSBML}(),
                           Dict{String, FunctionSBML}(), # SBML reactions
                           Dict{String, String}(), # Algebraic rules
                           Dict{String, String}(), # Generated id:s
                           Dict{String, String}(), # Piecewise to ifelse_expressions
                           Dict{String, String}(), # Ifelse to bool expression
                           Vector{String}(undef, 0), # Rate rule variables
                           Vector{String}(undef, 0), # Assignment rule variables
                           Vector{String}(undef, 0), # Algebraic rule variables
                           Vector{String}(undef, 0), # Species_appearing in reactions
                           Vector{String}(undef, 0),
                           conversion_factor,
                           specie_reference_ids,
                           Vector{String}(undef, 0),
                           libsbml_rule_variables,
                           reactionids) # Variables with piecewise
    return model_SBML
end
function ModelSBML(libsbml_model::SBML.Model)::ModelSBML
    conversion_factor = _parse_variable(libsbml_model.conversion_factor)
    name = _parse_variable(libsbml_model.name; default="SBML_model")
    name = replace(name, " " => "_")

    # Specie reference ids can sometimes appear in math expressions, where they should
    # be replaced. Precomputing the ids save computational time. Similar hold for rule
    # variables
    specie_reference_ids = get_specie_reference_ids(libsbml_model)
    libsbml_rule_variables = _get_sbml_rules_variables(libsbml_model)
    if isempty(libsbml_model.reactions)
        reactionids = String[]
    else
        reactionids = collect(keys(libsbml_model.reactions))
    end
    return ModelSBML(name, specie_reference_ids=specie_reference_ids, conversion_factor=conversion_factor, libsbml_rule_variables=libsbml_rule_variables, reactionids=reactionids)
end

struct ModelSBMLSystem
    species::String
    specie_map::String
    variables::String
    parameters::String
    parameter_map::String
    reactions::String
    has_species::Bool
    all_integer_S::Bool
end

struct SBMLSupport <: Exception
    var::String
end

mutable struct MathSBML
    formula::String
    math_idents::Vector{String}
    tmp_arg::String
    args::Vector{String}
    fns::Vector{String}
    user_fns::Vector{String}
    has_reaction_ids::Bool
    has_rateOf::Bool
end
