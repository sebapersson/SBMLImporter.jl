function parse_SBML_reactions!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    for (id, reaction) in libsbml_model.reactions
        # Process kinetic math into Julia syntax
        _formula, math_idents = parse_SBML_math(reaction.kinetic_math)

        # Check if reactionId or assignment-rule variable appear in the kinetic math,
        # if they do they are a part of the math_idents, and need to be replaced
        # in later processing steps
        has_assignment_rule_variable::Bool = false
        has_reaction_id::Bool = false
        for math_ident in math_idents
            if haskey(libsbml_model.reactions, math_ident)
                has_reaction_id = true
            end
            if math_ident ∈ model_SBML.assignment_rule_variables
                has_assignment_rule_variable = true
            end
        end

        # Add values for potential kinetic parameters (where-statements)
        for (parameter_id, parameter) in reaction.kinetic_parameters
            _formula = replace_variable(_formula, parameter_id, string(parameter.value))
        end

        # Capture potential piecewise
        if occursin("piecewise(", _formula)
            _formula = piecewise_to_ifelse(_formula, model_SBML, libsbml_model)
        end

        # Replace SBML functions, rescale species properly etc...
        formula = process_SBML_str_formula(_formula, model_SBML, libsbml_model,
                                           check_scaling = true)
        reactants = Vector{String}(undef, length(reaction.reactants))
        reactants_stoichiometry = similar(reactants)
        products = Vector{String}(undef, length(reaction.products))
        products_stoichiometry = similar(products)

        stoichiometry_mass_action::Bool = true
        for (i, reactant) in pairs(reaction.reactants)
            push!(model_SBML.species_in_reactions, reactant.species)
            if model_SBML.species[reactant.species].boundary_condition == true
                reactants_stoichiometry[i] = "nothing"
                reactants[i] = "nothing"
                continue
            end
            stoichiometry, _stoichiometry_mass_action = parse_stoichiometry(reactant,
                                                                            model_SBML)
            reactants[i] = reactant.species
            stoichiometry_mass_action = stoichiometry_mass_action == true &&
                                        _stoichiometry_mass_action == false ? false :
                                        stoichiometry_mass_action
            compartment_scaling = get_compartment_scaling(reactant.species, formula,
                                                          model_SBML)
            model_SBML.species[reactant.species].formula *= " - " * stoichiometry *
                                                            compartment_scaling * "(" *
                                                            formula * ")"

            # In case a specie is given in unit concentration we must scale the stoichiometry
            reactants_stoichiometry[i] = stoichiometry *
                                         get_compartment_scaling(reactant.species, formula,
                                                                 model_SBML;
                                                                 stoichiometry = true)
        end

        for (i, product) in pairs(reaction.products)
            push!(model_SBML.species_in_reactions, product.species)
            if model_SBML.species[product.species].boundary_condition == true
                products_stoichiometry[i] = "nothing"
                products[i] = "nothing"
                continue
            end
            stoichiometry, _stoichiometry_mass_action = parse_stoichiometry(product,
                                                                            model_SBML)
            stoichiometry_mass_action = stoichiometry_mass_action == true &&
                                        _stoichiometry_mass_action == false ? false :
                                        stoichiometry_mass_action
            compartment_scaling = get_compartment_scaling(product.species, formula,
                                                          model_SBML)
            model_SBML.species[product.species].formula *= " + " * stoichiometry *
                                                           compartment_scaling * "(" *
                                                           formula * ")"

            products[i] = product.species
            # See comment above regarding scaling for species given in concentration
            products_stoichiometry[i] = stoichiometry *
                                        get_compartment_scaling(product.species, formula,
                                                                model_SBML;
                                                                stoichiometry = true)
        end

        model_SBML.reactions[id] = ReactionSBML(id, formula, products,
                                                products_stoichiometry, reactants,
                                                reactants_stoichiometry,
                                                stoichiometry_mass_action,
                                                has_assignment_rule_variable,
                                                has_reaction_id)
    end

    # Species given via assignment rules, or initial assignments which only affect stoichiometry
    # are species that should not be included down the line in the model, hence they are
    # here removed from the model
    remove_stoichiometry_math_from_species!(model_SBML, libsbml_model)
end

function get_compartment_scaling(specie::String, formula::String, model_SBML::ModelSBML;
                                 stoichiometry::Bool = false)::String

    # The case of the specie likelly being given via an algebraic rule
    if isempty(formula)
        return stoichiometry == false ? "*" : ""
    end

    if model_SBML.species[specie].only_substance_units == true
        return stoichiometry == false ? "*" : ""
    end

    if model_SBML.species[specie].unit == :Amount
        return stoichiometry == false ? "*" : ""
    end

    if model_SBML.species[specie].unit == :Concentration
        if stoichiometry == false
            return "/" * model_SBML.species[specie].compartment * "*"
        else
            return "/" * model_SBML.species[specie].compartment
        end
    end
end

function parse_stoichiometry(specie_reference::SBML.SpeciesReference,
                             model_SBML::ModelSBML)::Tuple{String, Bool}
    mass_action_stoichiometry = isnothing(specie_reference.id)
    if !isnothing(specie_reference.id)
        if haskey(model_SBML.generated_ids, specie_reference.id)
            stoichiometry = model_SBML.generated_ids[specie_reference.id]
            return stoichiometry, mass_action_stoichiometry

            # Two following special cases where the stoichiometry is given by another variable in the model
        elseif specie_reference.id ∈ model_SBML.rate_rule_variables
            return specie_reference.id, mass_action_stoichiometry

        elseif !isempty(model_SBML.events) && any(occursin.(specie_reference.id,
                             reduce(vcat, [e.formulas for e in values(model_SBML.events)])))
            return specie_reference.id, mass_action_stoichiometry
        end

        stoichiometry = specie_reference.id
        # Here the stoichiometry is given as an assignment rule which has been interpreted as an additional model specie,
        # so the value is taken from the initial value
        if haskey(model_SBML.species, stoichiometry)
            stoichiometry = model_SBML.species[stoichiometry].initial_value
            if is_number(stoichiometry)
                stoichiometry = isnothing(stoichiometry) || stoichiometry == "nothing" ?
                                "1" : stoichiometry
            end
            # Can be nested 1 level
            if haskey(model_SBML.species, stoichiometry) &&
               model_SBML.species[stoichiometry].assignment_rule == true
                stoichiometry = model_SBML.species[stoichiometry].initial_value
            end

            return stoichiometry, mass_action_stoichiometry
        end

        # Last case where stoichiometry is not referenced anywhere in the model assignments, rules etc..., assign
        # to default value
        stoichiometry = string(specie_reference.stoichiometry)

    else
        stoichiometry = isnothing(specie_reference.stoichiometry) ? "1" :
                        string(specie_reference.stoichiometry)
        stoichiometry = stoichiometry[1] == '-' ? "(" * stoichiometry * ")" : stoichiometry
    end

    _stoichiometry = isnothing(stoichiometry) || stoichiometry == "nothing" ? "1" :
                     stoichiometry
    return _stoichiometry, mass_action_stoichiometry
end

function remove_stoichiometry_math_from_species!(model_SBML::ModelSBML,
                                                 libsbml_model::SBML.Model)::Nothing
    for (id, reaction) in libsbml_model.reactions
        specie_references = vcat([reactant for reactant in reaction.reactants],
                                 [product for product in reaction.products])
        for specie_reference in specie_references
            if haskey(model_SBML.generated_ids, specie_reference.id) ||
               isnothing(specie_reference.id)
                continue
            end

            if specie_reference.id ∈ model_SBML.rate_rule_variables
                if haskey(libsbml_model.initial_assignments, specie_reference.id)
                    continue
                end
                stoichiometry_t0 = isnothing(specie_reference.stoichiometry) ? "1.0" :
                                   string(specie_reference.stoichiometry)
                model_SBML.species[specie_reference.id].initial_value = stoichiometry_t0
                continue
            end

            if !isempty(model_SBML.events) && any(occursin.(specie_reference.id,
                             reduce(vcat, [e.formulas for e in values(model_SBML.events)])))
                continue
            end

            # An artifact from how the stoichiometry is procssed as assignment rule
            # or initial assignment
            if haskey(model_SBML.species, specie_reference.id)
                delete!(model_SBML.species, specie_reference.id)
            end
            if specie_reference.id ∈ model_SBML.assignment_rule_variables
                filter!(x -> x != specie_reference.id, model_SBML.assignment_rule_variables)
            end
        end
    end
    return nothing
end
