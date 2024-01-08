"""
    SBML_to_ReactionSystem(path_SBML::AbstractString;
                           ifelse_to_callback::Bool=true,
                           write_to_file::Bool=false, 
                           verbose::Bool=true, 
                           return_all::Bool=false, 
                           model_as_string::Bool=false)

Parse an SBML model into a Catalyst `ReactionSystem` and potentially convert events/piecewise to callbacks.

For information on simulating the `ReactionSystem`, refer to the documentation.

For converting the SBML model directly into a ModelingToolkit `ODESystem` see the function `SBML_to_ReactionSystem`.

For testing `path_SBML` can be the model as a string if `model_as_string=true`.

!!! note
    The number of returned arguments depends on whether the SBML model has events and/or piecewise expressions (see below).

## Arguments
- `path_SBML`: File path to a valid SBML file (level 2 or higher).
- `ifelse_to_callback=true`: Whether to rewrite `ifelse` (piecewise) expressions to callbacks; recommended for performance.
- `write_to_file=false`: Whether to write the parsed SBML model to a Julia file in the same directory as the SBML file.
- `verbose=true`: Whether or not to display information on the number of return arguments.
- `return_all=true`: Whether or not to return all possible arguments (see below), regardless of whether the model has events.
- `model_as_string=false` : Whether or not the model (`path_SBML`) is provided as a string, for testing.

## Returns
- `rn`: A Catalyst `ReactionSystem` that for example can be converted into an `ODEProblem` and solved.
- `specie_map`: A species map setting initial values; together with the `ReactionSystem`, it can be converted into an `ODEProblem`.
- `parameter_map` A parameter map setting parameter values; together with the `ReactionSystem`, it can be converted into an `ODEProblem`.
- `cbset` - **only for models with events/piecewise expressions**: Callbackset (events) for the model.
- `get_tstops`- **Only for models with events/piecewise expressions**: Function computing time stops for discrete callbacks in the `cbset`.

## Examples
```julia
# Import and simulate model without events
using SBMLImporter
rn, specie_map, parameter_map = SBML_to_ReactionSystem(path_SBML)
sys = convert(ODESystem, rn)

using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Solve ODE with Rodas5P solver
sol = solve(prob, Rodas5P())
```
```julia
# Import a model with events
using SBMLImporter
rn, specie_map, parameter_map, cb, get_tstops = SBML_to_ReactionSystem(path_SBML)
sys = convert(ODESystem, rn)

using OrdinaryDiffEq
tspan = (0, 10.0)
prob = ODEProblem(sys, specie_map, tspan, parameter_map, jac=true)
# Compute event times
tstops = get_tstops(prob.u0, prob.p)
sol = solve(prob, Rodas5P(), tstops=tstops, callback=callbacks)
```
"""                 
function SBML_to_ReactionSystem(path_SBML::T;
                                ifelse_to_callback::Bool=true,
                                inline_assignment_rules::Bool=true,
                                write_to_file::Bool=false, 
                                verbose::Bool=true, 
                                ret_all::Bool=false, 
                                model_as_string::Bool=false) where T <: AbstractString

    # Intermediate model representation of a SBML model which can be processed into
    # an ODESystem
    model_SBML = build_SBML_model(path_SBML; ifelse_to_callback=ifelse_to_callback, model_as_string=model_as_string, 
                                  inline_assignment_rules=inline_assignment_rules)
    rule_variables = unique(vcat(model_SBML.rate_rule_variables, model_SBML.assignment_rule_variables,
                                 model_SBML.algebraic_rule_variables))

    # If model is written to file save it in the same directory as the SBML-file
    if model_as_string == false
        dir_save = joinpath(splitdir(path_SBML)[1], "SBML")
    else
        dir_save = joinpath(@__DIR__, "SBML")
    end
    if write_to_file == true && !isdir(dir_save)
        mkdir(dir_save)
    end
    path_save_model = joinpath(dir_save, model_SBML.name * ".jl")

    # Build the ReactionSystem. Must be done via Meta-parse, because if a function is used 
    # via @RuntimeGeneratedFunction runtime is very slow for large models
    (_species_write, _specie_map_write, _variables_write, 
     _parameters_write, _parameter_map_write, _reactions_write, 
     no_species, integer_stoichiometries) = _reactionsystem_from_SBML(model_SBML, rule_variables)
    # The model can have have only species or only variables or both. If it has variables they 
    # are given via SBML rules
    eval(Meta.parse("ModelingToolkit.@variables t"))
    eval(Meta.parse("D = Differential(t)"))
    sps = no_species ? Any[] : eval(Meta.parse(split(_species_write, "\n")[2]))
    vs = isempty(rule_variables) ? Any[] : eval(Meta.parse(_variables_write))
    if isempty(rule_variables)
        sps_arg = sps
    elseif no_species == false
        sps_arg = [sps; vs]
    else
        sps_arg = vs
    end                                  
    # Parameters can not be an empty collection
    if _parameters_write != "\tps = Catalyst.@parameters "
        ps = eval(Meta.parse(_parameters_write))
    else
        ps = Any[]
    end
    _reactions = eval(Meta.parse(_reactions_write))
    combinatoric_ratelaws = integer_stoichiometries ? true : false
    # Build reaction system from its components
    reaction_system = Catalyst.ReactionSystem(_reactions, t, sps_arg, ps; name=Symbol(model_SBML.name), combinatoric_ratelaws=combinatoric_ratelaws)
    specie_map = eval(Meta.parse(_specie_map_write))
    parameter_map = eval(Meta.parse(_parameter_map_write))

    # Build callback functions 
    cbset, callback_str = create_callbacks_SBML(reaction_system, model_SBML, model_SBML.name)

    # if model is written to file write the callback
    if write_to_file == true
        path_save = joinpath(dir_save, model_SBML.name * "_callbacks.jl")
        io = open(path_save, "w")
        write(io, callback_str)
        close(io)
        _ = reactionsystem_to_string(_species_write, _specie_map_write, _variables_write, 
                                     _parameters_write, _parameter_map_write, _reactions_write, 
                                     no_species, integer_stoichiometries, write_to_file, path_save_model, 
                                     rule_variables, model_SBML)
    end

    if ret_all == true
        return reaction_system, specie_map, parameter_map, cbset
    end

    # Depending on model return what is needed to perform forward simulations
    if isempty(model_SBML.events) && isempty(model_SBML.ifelse_bool_expressions)
        return reaction_system, specie_map, parameter_map
    end
    
    verbose && @info "SBML model with events - output returned as odesys, specie_map, parameter_map, cbset\nFor how to simulate model see documentation"
    return reaction_system, specie_map, parameter_map, cbset
end



function _reactionsystem_from_SBML(model_SBML::ModelSBML, 
                                   rule_variables::Vector{String})::Tuple{String, String, 
                                                                          String, String, 
                                                                          String, String, 
                                                                          Bool, Bool}

    # Check if model is empty of derivatives if the case add dummy state to be able to
    # simulate the model
    if ((isempty(model_SBML.species) || sum([!s.assignment_rule for s in values(model_SBML.species)]) == 0) &&
        (isempty(model_SBML.parameters) || sum([p.rate_rule for p in values(model_SBML.parameters)]) == 0) &&
        (isempty(model_SBML.compartments) || sum([c.rate_rule for c in values(model_SBML.compartments)]) == 0))

        model_SBML.species["foo"] = SpecieSBML("foo", false, false, "1.0", "0.0", "1.0", "", :Amount,
                                                  false, false, false, false)
    end

    # Setup Catalyst ReactionNetwork
    _species_write, _species_write_array, _specie_map_write = SBMLImporter.get_specie_map(model_SBML, reaction_system=true)
    _parameters_write, _parameters_write_array, _parameter_map_write = SBMLImporter.get_parameter_map(model_SBML, reaction_system=true)

    # In case a specie (or parameter) appear as a rate-rule, algebraic or assignment rule they need to be treated as
    # MTK variable for the the downstream processing. This might turn the species block empty, then it must be removed
    _variables_write = "\tvs = ModelingToolkit.@variables"
    filter!(x -> x ∉ keys(model_SBML.generated_ids), rule_variables)
    for variable in rule_variables
        _species_write = replace(_species_write, " " * variable * "(t)" => "")
        _variables_write *= " " * variable * "(t)"
    end
    if !isempty(rule_variables)
        no_species = all([x ∈ rule_variables for x in keys(model_SBML.species)])
    else
        no_species = false
    end

    # Reaction stoichiometry and propensities
    integer_stoichiometries::Bool = true
    _reactions_write = "\t_reactions = [\n"
    for (id, r) in model_SBML.reactions

        # Can happen for models with species that are boundary conditions, such species 
        # do not take part in reactions 
        if all(r.products .== "nothing") && all(r.reactants .== "nothing")
            continue
        end

        reactants_stoichiometries, reactants, integer_stoichiometries1 = get_reaction_side(r, :Reactants, model_SBML)
        products_stoichiometries, products, integer_stoichiometries2 = get_reaction_side(r, :Products, model_SBML)
        propensity = r.kinetic_math
        if reaction_is_mass_action(r, model_SBML) == true && integer_stoichiometries1 && integer_stoichiometries2
            _reactions_write *= ("\t\tSBMLImporter.update_rate_reaction(Catalyst.Reaction(" * propensity * ", " *
                                reactants * ", " * products * ", " *
                                reactants_stoichiometries * ", " * products_stoichiometries * "; only_use_rate=false)),\n")
        else
            _reactions_write *= ("\t\tCatalyst.Reaction(" * propensity * ", " *
                                reactants * ", " * products * ", " *
                                reactants_stoichiometries * ", " * products_stoichiometries * "; only_use_rate=true),\n")
        end
                        
        # If it has already been assigned false we know that all stoichiometries are not
        # integer numbers, and if either of integer_stoichiometries are false integer_stoichiometries
        # should become false
        if integer_stoichiometries == true
            integer_stoichiometries = !any([integer_stoichiometries1, integer_stoichiometries2] .== false)
        end
    end

    # Rules are directly encoded into the Catalyst.Reaction vector
    for variable in unique(vcat(model_SBML.rate_rule_variables, model_SBML.assignment_rule_variables))
        if variable ∈ keys(model_SBML.species)
            @unpack formula, assignment_rule, rate_rule = model_SBML.species[variable]
        elseif variable ∈ keys(model_SBML.parameters)
            @unpack formula, assignment_rule, rate_rule = model_SBML.parameters[variable]
        elseif variable ∈ keys(model_SBML.compartments)
            @unpack formula, assignment_rule, rate_rule = model_SBML.compartments[variable]
        else
            continue
        end
        if rate_rule == true
            _reactions_write *= "\t\tD(" * variable * ") ~ " * formula * ",\n"
        elseif assignment_rule == true
            _reactions_write *= "\t\t" * variable * " ~ " * formula * ",\n"
        end
    end
    for formula in values(model_SBML.algebraic_rules)
        _reactions_write *= "\t\t" * formula * ",\n"
    end
    _reactions_write *= "\t]\n"

    return (_species_write, _specie_map_write, _variables_write, 
            _parameters_write, _parameter_map_write, _reactions_write, 
            no_species, integer_stoichiometries)

end


function reactionsystem_to_string(_species_write::String, _specie_map_write::String, 
                                  _variables_write::String, _parameters_write::String, 
                                  _parameter_map_write::String, _reactions_write::String, 
                                  no_species::Bool, integer_stoichiometries::Bool, 
                                  write_to_file::Bool, path_save_model::String, 
                                  rule_variables::Vector{String}, model_SBML::ModelSBML)::String

    # ReactionSystem
    combinatoric_ratelaws_arg = integer_stoichiometries ? "true" : "false"
    if isempty(rule_variables)
        sps_arg = "sps"
    elseif no_species == false
        sps_arg = "[sps; vs]"
    else
        sps_arg = "vs"
        _species_write = replace(_species_write, "sps = Catalyst.@species" => "")
    end                                  
    # Parameters might be an empty set
    if _parameters_write != "\tps = Catalyst.@parameters "
        _rn_write = "\trn = Catalyst.ReactionSystem(reactions, t, $sps_arg, ps; name=Symbol(\"" * model_SBML.name * "\"), combinatoric_ratelaws=$combinatoric_ratelaws_arg)"
    else
        _rn_write = "\trn = Catalyst.ReactionSystem(reactions, t, $sps_arg, Any[]; name=Symbol(\"" * model_SBML.name * "\"), combinatoric_ratelaws=$combinatoric_ratelaws_arg)"
    end

    # Create a function returning the ReactionSystem, specie-map, and parameter-map
    _function_write = "function get_reaction_system(foo)\n\n"
    _function_write *= _species_write * "\n"
    if _variables_write != "\tvs = ModelingToolkit.@variables"
        _function_write *= _variables_write * "\n"
    end
    if _parameters_write != "\tps = Catalyst.@parameters "
        _function_write *= _parameters_write * "\n\n"
    end
    _function_write *= "\tD = Differential(t)\n\n"
    _function_write *= _reactions_write * "\n\n"
    _function_write *= _rn_write * "\n\n"
    _function_write *= _specie_map_write * "\n"
    _function_write *= _parameter_map_write * "\n"
    _function_write *= "\treturn rn, specie_map, parameter_map\nend"

    # In case user request file to be written
    if write_to_file == true
        open(path_save_model, "w") do f
            write(f, _function_write)
        end
    end

    return _function_write
end



function get_reaction_side(r::ReactionSBML, which_side::Symbol, 
                           model_SBML::ModelSBML)::Tuple{String, String, Bool}

    if which_side === :Reactants
        species, stoichiometries = r.reactants, r.reactants_stoichiometry
    elseif which_side === :Products
        species, stoichiometries = r.products, r.products_stoichiometry
    end

    # Case where we go from ϕ -> prod (or reverse)
    if isempty(species)
        return "nothing", "nothing", true
    end
    # Edge case for boundary condition and rate rule
    if length(species) == 1 && species[1] == "nothing"
        return "nothing", "nothing", true
    end

    # Process the _toichiometry for the reaction, species vector is not 
    # required to be unique hence potential double counting must be 
    # considered 
    k = 1
    _stoichiometries = Vector{String}(undef, length(filter(x -> x != "nothing", unique(species))))
    _species_parsed = fill("", length(filter(x -> x != "nothing", unique(species))))

    # Happens when all reactants or products are boundary condition
    if isempty(_species_parsed)
        return "nothing", "nothing", true
    end

    integer_stoichiometry::Bool = true
    for i in eachindex(species)

        # Happens when a specie is a boundary condition, it should not be involved in the reaction and 
        # should affect reaction dynamics
        if species[i] == "nothing"
            continue
        end

        stoichiometry, _integer_stoichiometry = parse_stoichiometry_reaction_system(stoichiometries[i])

        # SBML models can have conversion factors that scale stoichiometry
        if isempty(model_SBML.species[species[i]].conversion_factor) && isempty(model_SBML.conversion_factor)
            _stoichiometry = stoichiometry 
            if integer_stoichiometry == true
                integer_stoichiometry = _integer_stoichiometry
            end
        else
            cv_specie = model_SBML.species[species[i]].conversion_factor
            cv = isempty(cv_specie) ? model_SBML.conversion_factor : cv_specie
            _stoichiometry = stoichiometry * "*" * cv 
            integer_stoichiometry = false
        end

        # Species are allowed to be repeated in SBML reactions, this is not 
        # allowed in Catalyst, therefore stoichiometry for these repititions 
        # are added up 
        if species[i] ∈ _species_parsed
            _i = findfirst(x -> x == species[i], _species_parsed)
            _stoichiometries[_i] *= "+" * _stoichiometry
        else
            _species_parsed[k] = species[i]
            _stoichiometries[k] = _stoichiometry
            k += 1
        end

    end
    _stoichiometries_str = "[" * prod([s * ", " for s in _stoichiometries])[1:end-2] * "]"
    _species_str = "[" * prod([s * ", " for s in _species_parsed])[1:end-2] * "]"

    return _stoichiometries_str, _species_str, integer_stoichiometry
end


# If possible parse stoichiometry to an integer
function parse_stoichiometry_reaction_system(stoichiometry::String)::Tuple{String, Bool}
    if !isnothing(tryparse(Float64, stoichiometry))
        _stoichiometry = parse(Float64, stoichiometry)
        try
            return string(Int64(_stoichiometry)), true
        catch
            return string(_stoichiometry), false
        end
    else
        return stoichiometry, false
    end
end


"""
    function reaction_is_mass_action(r::ReactionSBML, model_SBML::ModelSBML)::Bool

Check if a Catalyst recation should be converted to a mass-action reaction, occurs if:

* Involved species has only_substance_unit=true
* Propensity does not include time t, or depend on rate-rule or assignment rule variables
* Stoichiometry is mass-action. In case reactants or products have a specie-id then the reaction 
  does not have to be mass-action following SBML standard
"""
function reaction_is_mass_action(r::ReactionSBML, model_SBML::ModelSBML)::Bool

    if r.stoichiometry_mass_action == false
        return false
    end

    for reactant in r.reactants
        if reactant == "nothing"
            continue
        end
            
        if model_SBML.species[reactant].only_substance_units == false
            return false
        end
    end
    for product in r.products
        if product == "nothing"
            continue
        end

        if model_SBML.species[product].only_substance_units == false
            return false
        end
    end

    # Check that no rule variables appear in formula
    formula = r.kinetic_math
    rule_variables = unique(vcat(model_SBML.rate_rule_variables, model_SBML.assignment_rule_variables,
                                 model_SBML.algebraic_rule_variables))
    for rule_variable in rule_variables
        if SBMLImporter.replace_variable(formula, rule_variable, "") != formula
            return false
        end
    end

    return true
end


function update_rate_reaction(rx; combinatoric_ratelaw::Bool=true)
    @set rx.rate = (rx.rate^2) / Catalyst.oderatelaw(rx; combinatoric_ratelaw=combinatoric_ratelaw)
end