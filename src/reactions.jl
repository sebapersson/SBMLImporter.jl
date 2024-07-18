function parse_reactions!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (id, reaction) in libsbml_model.reactions
        propensity, math_expression = _parse_reaction_formula(reaction, model_SBML,
                                                              libsbml_model)

        reactants = _get_reaction_species(reaction, model_SBML, :reactants)
        reactants_cs = _get_compartment_scalings(reactants, propensity, model_SBML)
        reactants_s, reactants_massaction = _get_stoichiometries(reaction, reactants,
                                                                 model_SBML, :reactants)
        for (i, reactant) in pairs(reactants)
            reactant == "nothing" && continue
            _update_ode!(model_SBML.species[reactant], reactants_s[i], reactants_cs[i],
                         propensity, :reactant)
            reactants_s[i] = _template_stoichiometry(reactants_s[i], reactants_cs[i])
            _add_ident_info!(model_SBML.species[reactant], math_expression, model_SBML)
        end

        products = _get_reaction_species(reaction, model_SBML, :products)
        products_cs = _get_compartment_scalings(products, propensity, model_SBML)
        products_s, products_massaction = _get_stoichiometries(reaction, products,
                                                               model_SBML, :products)
        for (i, product) in pairs(products)
            product == "nothing" && continue
            _update_ode!(model_SBML.species[product], products_s[i], products_cs[i],
                         propensity, :product)
            products_s[i] = _template_stoichiometry(products_s[i], products_cs[i])
            _add_ident_info!(model_SBML.species[product], math_expression, model_SBML)
        end

        # Storing whether kinetic_math has assignment rules and reaction-ids is important
        # for faster downstream processes, as these can be replaced
        assignment_rules = _has_assignment_rule_ident(math_expression.math_idents,
                                                      model_SBML)
        have_ridents = math_expression.has_reaction_ids
        have_rateOf = "rateOf" in math_expression.fns
        have_specieref = _has_specieref(math_expression.math_idents, model_SBML)
        stoichiometry_massaction = all([reactants_massaction, products_massaction])
        model_SBML.reactions[id] = ReactionSBML(id, propensity, products, products_s,
                                                reactants, reactants_s,
                                                stoichiometry_massaction, assignment_rules,
                                                have_ridents, have_rateOf, have_specieref)
        for specie in Iterators.flatten((reactants, products))
            push!(model_SBML.species_in_reactions, specie)
        end
    end

    # If a ModelSBML specie appears in a SpeciesReference.id this is because of a
    # StoichiometryMath. In libsbml parsing the StoichiometryMath is promoted to an
    # assignment rule. As this rule variable does not correspond to any model variable,
    # followning the standard a new specie is created. However, this is just an artifact,
    # and the user never meant to create a specie. Hence, this step removes the
    # StoichiometryMath species from the model
    _remove_stoichiometry_math!(model_SBML, libsbml_model)
    return nothing
end

function _parse_reaction_formula(reaction::SBML.Reaction, model_SBML::ModelSBML,
                                 libsbml_model::SBML.Model)::Tuple{String, MathSBML}
    math_expression = parse_math(reaction.kinetic_math, libsbml_model)
    formula = math_expression.formula
    # SBML where statements, that can occur in reaction
    for (parameter_id, parameter) in reaction.kinetic_parameters
        formula = _replace_variable(formula, parameter_id, string(parameter.value))
    end

    math_expression.formula = formula
    formula = _process_formula(math_expression, model_SBML)
    return formula, math_expression
end

function _get_reaction_species(reaction::SBML.Reaction, model_SBML::ModelSBML,
                               which_side::Symbol)::Vector{String}
    species = which_side == :reactants ? reaction.reactants : reaction.products

    reaction_species = fill("", length(species))
    for (i, specie) in pairs(species)
        id = specie.species
        if model_SBML.species[id].boundary_condition == true
            reaction_species[i] = "nothing"
            continue
        end
        reaction_species[i] = id
    end
    return reaction_species
end

function _get_stoichiometries(reaction::SBML.Reaction, species_id::Vector{String},
                              model_SBML::ModelSBML,
                              which_side::Symbol)::Tuple{Vector{String}, Bool}
    species = which_side == :reactants ? reaction.reactants : reaction.products

    massaction::Bool = true
    stoichiometries = similar(species_id)
    for (i, specie) in pairs(species)
        if species_id[i] == "nothing"
            stoichiometries[i] = "nothing"
            continue
        end
        stoichiometries[i], massaction = _parse_stoichiometry(specie, model_SBML)
    end
    return stoichiometries, massaction
end

function _get_compartment_scalings(species::Vector{String}, propensity::String,
                                   model_SBML::ModelSBML)::Vector{String}
    compartment_scalings = similar(species)
    for (i, specie) in pairs(species)
        if isempty(propensity)
            # Likelly specie is given by algebraic rule
            compartment_scalings[i] = "1.0"
            continue
        end
        if specie == "nothing"
            compartment_scalings[i] = "1.0"
            continue
        end
        if model_SBML.species[specie].only_substance_units == true
            compartment_scalings[i] = "1.0"
            continue
        end
        if model_SBML.species[specie].unit == :Amount
            compartment_scalings[i] = "1.0"
            continue
        end
        @assert model_SBML.species[specie].unit==:Concentration "Problem parsing compartment for reactions"
        compartment_scalings[i] = _apply(/, "1.0", model_SBML.species[specie].compartment)
    end
    return compartment_scalings
end

function _parse_stoichiometry(specie::SBML.SpeciesReference,
                              model_SBML::ModelSBML)::Tuple{String, Bool}
    # specie.id (that can be a SBML variable) has priority over specie.stoichiometry.
    if isnothing(specie.id)
        s = _parse_variable(specie.stoichiometry; default = "1")
        s = s[1] == '-' ? "(" * s * ")" : s
        return s, true
    end

    # Stoichiometry can be given by other model variables that already have a formula,
    # such as rate rules.
    massaction = false
    if haskey(model_SBML.generated_ids, specie.id)
        s = model_SBML.generated_ids[specie.id]
        return s, massaction
    end
    if specie.id in model_SBML.rate_rule_variables
        return specie.id, massaction
    end
    if _is_event_assigned(specie.id, model_SBML) == true
        return specie.id, massaction
    end

    # Here the stoichiometry is given as an assignment rule, that has been added as a
    # model specie during rule parsing
    if haskey(model_SBML.species, specie.id)
        s = model_SBML.species[specie.id].initial_value
        # 1 level of nesting allowed
        if haskey(model_SBML.species, s) && model_SBML.species[s].assignment_rule == true
            s = model_SBML.species[s].initial_value
        end
        return s, massaction
    end

    # At this point the model has a SpeciesReference.id that does not correspond to any
    # model variables. For example happens with rules that have no math-child, in this case
    # the specie.stoichiometry should be used
    s = _parse_variable(specie.stoichiometry; default = "1")
    return s, massaction
end

function _update_ode!(specie::SpecieSBML, s::String, c_scaling::String, propensity::String,
                      which_side::Symbol)::Nothing
    ode_formula = _template_ode_reaction(s, c_scaling, propensity, which_side)
    if isempty(specie.formula)
        specie.formula = ode_formula
    else
        specie.formula = _apply(+, specie.formula, ode_formula)
    end
    return nothing
end

function _remove_stoichiometry_math!(model_SBML::ModelSBML,
                                     libsbml_model::SBML.Model)::Nothing
    for reaction in values(libsbml_model.reactions)
        species_ref = vcat([r for r in reaction.reactants], [p for p in reaction.products])
        for specie_ref in species_ref
            if isnothing(specie_ref.id)
                continue
            end
            if haskey(model_SBML.generated_ids, specie_ref.id)
                continue
            end
            if specie_ref.id in model_SBML.rate_rule_variables
                continue
            end
            if _is_event_assigned(specie_ref.id, model_SBML)
                continue
            end

            if haskey(model_SBML.species, specie_ref.id)
                delete!(model_SBML.species, specie_ref.id)
            end
            if specie_ref.id in model_SBML.assignment_rule_variables
                filter!(x -> x != specie_ref.id, model_SBML.assignment_rule_variables)
            end
        end
    end
    return nothing
end
