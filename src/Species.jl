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
        # guaranteed zero.
        if boundary_condition == true || constant == true
            formula = "0.0"
        else
            formula = ""
        end
        model_SBML.species[specie_id] = SpecieSBML(specie_id, boundary_condition, constant, initial_value, formula, compartment, conversion_factor, unit, only_substance_units, false, false, false, false, false, false)
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

function adjust_for_dynamic_compartment!(model_SBML::ModelSBML)::Nothing
    # For species with boundary condition given by a rate-rule amount should stay constant
    # even if the compartment changes. It follows that amount n (amount), V (compartment)
    # and concentration (n/V) are related via:
    # dn/dt = d(n/V)/dt*V + n*(dV/dt)/V
    # Practically, this is handled by adding an additional model specie that describes the
    # concentration change, where the formula d(n/V)/dt is given by the rate-rule
    # divided by compartment
    for (specie_id, specie) in model_SBML.species
        compartment = model_SBML.compartments[specie.compartment]
        compartment.rate_rule == false && continue
        compartment.constant == true && continue
        specie.boundary_condition == false && continue
        specie.rate_rule == false && continue
        specie.only_substance_units == true && continue
        specie.unit != :Amount && continue

        conc_id = "__" * specie_id * "__conc__"
        conc_initial_value = specie.initial_value * "/" * compartment.name
        dcdt = specie.formula * "/" * compartment.name
        push!(model_SBML.rate_rule_variables, conc_id)
        model_SBML.species[conc_id] = SpecieSBML(conc_id, false, false, conc_initial_value, dcdt, compartment.name, specie.conversion_factor, :Concentration, false, false, true, false, specie.has_reaction_ids, specie.has_rateOf, specie.has_specieref)

        V, dVdt = compartment.name, compartment.formula
        specie.formula = dcdt * "*" * V * " + " * specie_id * "*" * dVdt * " / " * V
    end

    # When a specie is given in concentration but compartment changes, the concentration (c)
    # is given by:
    # dc/dt = (dn/dt - n*(dV/dt)/V)/V
    # Practically this is handled by adding an additional model specie that describes the
    # amount change, where the formula dn/dt is given by the Catalyst reactions, and the
    # conc. formula is given by the above ODE
    for (specie_id, specie) in model_SBML.species
        # To avoid that cases like above are parsed again
        occursin("__conc__", specie_id) && continue

        # Assignment rule can nest to rate-rules, in this case, the actual compartment
        # is the rate rule. If compartment.formula
        compartment = model_SBML.compartments[specie.compartment]
        if compartment.assignment_rule == true
            if compartment.formula in model_SBML.rate_rule_variables
                compartment = model_SBML.parameters[compartment.formula]
            elseif _contains_time_or_species(compartment.formula, model_SBML) == false
                continue
            end
        end
        compartment.constant == true && continue
        specie.unit != :Concentration && continue
        specie.rate_rule == true && continue
        specie.only_substance_units == true && continue
        # If a compartment is inconstant due to an event, it is (see events.jl) promoted
        # to rate-rule to not be simplified away by structurally_simplify. In this case
        # the formula = 0.0 and compartment changes are handled by the event parsing.
        if compartment.rate_rule == true
            compartment.formula == "0.0" && continue
        end

        V, dVdt = compartment.name, compartment.formula
        n_id = "__" * specie_id * "__amount__"
        n_initial_amount = specie.initial_value * "*" * V
        dndt = _get_amount_formula(specie, V)
        model_SBML.species[n_id] = SpecieSBML(n_id, false, false, n_initial_amount, dndt, V, specie.conversion_factor, :Amount, false, false, false, false, specie.has_reaction_ids, specie.has_rateOf, specie.has_specieref)

        # To enforce that dcdt is given by the ODE specie is promoted to rate-rule. Further
        # to obtain correct dn/dt formula, replace specie with n_id in reactions
        #V, dVdt = compartment.name, compartment.formula
        #dVdt = compartment.formula
        specie.formula = dndt * "/(" * V * ") - " * n_id * "/(" * V * ")^2*" * dVdt
        specie.rate_rule = true
        push!(model_SBML.rate_rule_variables, specie.name)
        for reaction in values(model_SBML.reactions)
            for (i, product) in pairs(reaction.products)
                product != specie.name && continue
                reaction.products[i] = n_id
                reaction.products_stoichiometry[i] *= "*" * V
            end
            for (i, reactant) in pairs(reaction.reactants)
                reactant != specie.name && continue
                reaction.reactants[i] = n_id
                reaction.reactants_stoichiometry[i] *= "*" * V
            end
        end
    end
    return nothing
end

function _get_amount_formula(specie::SpecieSBML, V::String)::String
    # If boundary condition is true; amount, not concentration should stay constant with
    # time
    if specie.boundary_condition == true
        return "0.0"
    end
    return isempty(specie.formula) ? "0.0" : "(" * specie.formula * ")" * V
end

function add_conversion_factor_ode!(model_SBML::ModelSBML, libsbml_model::SBML.Model)::Nothing
    # Conversion factor only apply to species changed via recations, not rules or
    # boundary conditions per SBML standard
    for (specie_id, specie) in model_SBML.species
        !haskey(libsbml_model.species, specie_id) && continue
        specie.assignment_rule == true && continue
        specie.rate_rule == true && continue
        specie.boundary_condition == true && continue
        isempty(specie.formula) && continue

        cf = _get_cf_scaling(specie, model_SBML)
        isempty(cf) && continue
        specie.formula = "(" * specie.formula * ") * " * cf
    end
    return nothing
end
