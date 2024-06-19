#=
    Functionality for parsing, and handling SBML species (e.g conversion factor etc...)
=#

function parse_species!(model_SBML::ModelSBML, libsbml_model::SBML.Model,
                        mass_action::Bool)::Nothing
    for (specie_id, specie) in libsbml_model.species
        if specie_id in FORBIDDEN_IDS
            throw(SBMLSupport("Specie name $(specie_id) is not allowed."))
        end
        @assert !isempty(specie.compartment) "Specie $specie_id does not have a compartment"

        # If mass_action=true the user enforces only_substance_units=true
        only_substance_units = _get_only_substance_units(specie, mass_action)
        unit = _get_unit(specie, only_substance_units)
        initial_value = _get_initial_value(specie, unit)
        conversion_factor = _parse_variable(specie.conversion_factor)
        constant = _parse_bool(specie.constant)
        boundary_condition = _parse_bool(specie.boundary_condition)
        compartment = specie.compartment
        # In case we have a constant or boundary condition the formula (derivative) is
        # guaranteed zero. TODO: I think this should be moved where formula parsing occurs
        if boundary_condition == true || constant == true
            formula = "0.0"
        else
            formula = ""
        end
        model_SBML.species[specie_id] = SpecieSBML(specie_id, boundary_condition, constant, initial_value, formula, compartment, conversion_factor, unit, only_substance_units, false, false, false, false)
    end
    return nothing
end

function _get_only_substance_units(specie::SBML.Species, mass_action::Bool)::Bool
    mass_action == true && return true
    return _parse_bool(specie.only_substance_units)
end

function _get_unit(specie::SBML.Species, only_substance_units::Bool)::Symbol
    @unpack initial_amount, initial_concentration, substance_units = specie
    # Per SBML standard if initial_concentration or initial_amount set the unit that is
    # the specie unit, otherwise, it is given by substance_units.
    if only_substance_units == true
        unit = :Amount
    elseif !isnothing(initial_concentration)
        unit = :Concentration
    elseif !isnothing(initial_amount)
        unit = :Amount
    else
        unit = substance_units == "substance" ? :Amount : :Concentration
    end
    return unit
end

function _get_initial_value(specie::SBML.Species, unit::Symbol)::String
    @unpack initial_concentration, only_substance_units, initial_amount = specie
    given_in_conc = !isnothing(initial_concentration)
    given_in_amount = !isnothing(initial_amount)
    if given_in_conc && unit == :Amount
        initial_value = string(initial_concentration) * "*" * specie.compartment
    elseif given_in_conc && unit == :Concentration
        initial_value = initial_concentration
    elseif given_in_amount && unit == :Concentration
        initial_value = string(initial_amount) * "/" * specie.compartment
    elseif given_in_amount && unit == :Amount
        initial_value = _parse_variable(initial_amount)
    # Defaults to zero
    else
        initial_value = "0.0"
    end
    return string(initial_value)
end

# TODO: This does not fit here. Dynamic compartments is its own thing
function adjust_for_dynamic_compartment!(model_SBML::ModelSBML)::Nothing

    #=
    The volume might change over time but the amount should stay constant, as we have a boundary condition
    for a specie given by a rate-rule. In this case it follows that amount n (amount), V (compartment) and conc.
    are related via the chain rule by:
    dn/dt = d(n/V)/dt*V + n*dV/dt/V
    =#
    for (specie_id, specie) in model_SBML.species

        # Specie with amount where amount should stay constant
        compartment = model_SBML.compartments[model_SBML.species[specie_id].compartment]
        if !(specie.unit == :Amount &&
             specie.rate_rule == true &&
             specie.boundary_condition == true &&
             specie.only_substance_units == false &&
             compartment.rate_rule == true)
            continue
        end

        if compartment.constant == true
            continue
        end

        # In this case must add additional variable for the specie concentration, to properly get the amount equation
        specie_conc_id = "__" * specie_id * "__conc__"
        initial_value_conc = model_SBML.species[specie_id].initial_value * "/" *
                             compartment.name
        formula_conc = model_SBML.species[specie_id].formula * "/" * compartment.name

        # Formula for amount specie. Treated as rate-rule as this is a feature we do not support
        # for mass-action models
        model_SBML.species[specie_id].formula = formula_conc * "*" * compartment.name *
                                                " + " * specie_id * "*" *
                                                compartment.formula * " / " *
                                                compartment.name

        # Add new conc. specie to model. See comment on rate-rule above
        model_SBML.species[specie_conc_id] = SpecieSBML(specie_conc_id, false, false,
                                                        initial_value_conc,
                                                        formula_conc, compartment.name,
                                                        specie.conversion_factor,
                                                        :Concentration, false, false, true,
                                                        false, specie.has_reaction_ids)
        push!(model_SBML.rate_rule_variables, specie_conc_id)
    end

    # When a specie is given in concentration, but the compartment concentration changes
    for (specie_id, specie) in model_SBML.species

        # To avoid that concentration species given as above are processed again
        if length(specie_id) ≥ 2 && specie_id[1:2] == "__"
            continue
        end

        compartment = model_SBML.compartments[model_SBML.species[specie_id].compartment]
        compartment_name = compartment.name
        if compartment.assignment_rule == true &&
           compartment.formula ∈ model_SBML.rate_rule_variables
            compartment = model_SBML.parameters[compartment.formula]
        end

        if !(specie.unit == :Concentration &&
             specie.only_substance_units == false &&
             compartment.constant == false)
            continue
        end
        # Rate rule has priority
        if specie_id ∈ model_SBML.rate_rule_variables
            continue
        end
        if !any(occursin.(keys(model_SBML.species), compartment.formula)) &&
           compartment.rate_rule == false
            continue
        end

        # Derivative and inital values newly introduced amount specie
        specie_amount_id = "__" * specie_id * "__amount__"
        initial_value_amount = specie.initial_value * "*" * compartment.name

        # If boundary condition is true only amount, not concentration should stay constant with
        # compartment size
        if specie.boundary_condition == true
            formula_amount = "0.0"
        else
            formula_amount = isempty(specie.formula) ? "0.0" :
                             "(" * specie.formula * ")" * compartment_name
        end

        # New formula for conc. specie
        specie.formula = formula_amount * "/(" * compartment_name * ") - " *
                         specie_amount_id * "/(" * compartment_name * ")^2*" *
                         compartment.formula
        specie.rate_rule = true
        push!(model_SBML.rate_rule_variables, specie.name)

        # Add new conc. specie to model
        model_SBML.species[specie_amount_id] = SpecieSBML(specie_amount_id, false, false,
                                                          initial_value_amount,
                                                          formula_amount, compartment_name,
                                                          specie.conversion_factor,
                                                          :Amount, false, false, false,
                                                          false, specie.has_reaction_ids)
        for (r_id, r) in model_SBML.reactions
            for i in eachindex(r.products)
                if r.products[i] == specie.name
                    r.products[i] = specie_amount_id
                    r.products_stoichiometry[i] *= "*" * compartment_name
                end
            end
            for i in eachindex(r.reactants)
                if r.reactants[i] == specie.name
                    r.reactants[i] = specie_amount_id
                    r.reactants_stoichiometry[i] *= "*" * compartment_name
                end
            end
        end
    end
    return nothing
end

# TODO : Should be in general equation parsing files
function adjust_conversion_factor!(model_SBML::ModelSBML,
                                   libsbml_model::SBML.Model)::Nothing
    for (specie_id, specie) in model_SBML.species
        if specie.assignment_rule == true
            continue
        end

        if !haskey(libsbml_model.species, specie_id)
            continue
        end

        # Conversion factors only affect species whose values are changed via reactions,
        # but not rate-rules
        if specie.rate_rule == true
            continue
        end

        # Boundary value species are not affected
        if specie.boundary_condition == true
            continue
        end

        # Zero change of rate for specie
        if isempty(specie.formula)
            continue
        end

        if !isempty(specie.conversion_factor)
            conversion_factor = specie.conversion_factor
        elseif !isnothing(libsbml_model.conversion_factor)
            conversion_factor = libsbml_model.conversion_factor
        else
            return nothing
        end

        specie.formula = "(" * specie.formula * ") * " * conversion_factor
    end

    return nothing
end
