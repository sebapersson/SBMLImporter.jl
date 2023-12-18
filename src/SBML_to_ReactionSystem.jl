function SBML_to_ReactionSystem(path_SBML::T;
                                ifelse_to_callback::Bool=true,
                                verbose::Bool=true,
                                ret_all::Bool=false,
                                model_as_string::Bool=false) where T <: AbstractString

    # Intermediate model representation of a SBML model which can be processed into
    # an ODESystem
    write_to_file = false
    model_SBML = build_SBML_model(path_SBML; ifelse_to_callback=ifelse_to_callback, model_as_string=model_as_string)

    # If model is written to file save it in the same directory as the SBML-file
    if model_as_string == false
        dir_save = joinpath(splitdir(path_SBML)[1], "SBML")
    else
        dir_save = joinpath(@__DIR__, "SBML")
    end
    if write_to_file == true && !isdir(dir_save)
        mkdir(dir_save)
    end

    # Build the ReactionSystem
    _model = reactionsystem_from_SBML(model_SBML)
    _get_reaction_system = @RuntimeGeneratedFunction(Meta.parse(_model))
    reaction_system, specie_map, parameter_map = _get_reaction_system("https://xkcd.com/303/") # Argument needed by @RuntimeGeneratedFunction

    # Build callback functions
    _callbacks, _tstops = create_callbacks_SBML(reaction_system, model_SBML, model_SBML.name)
    get_callbacks = @RuntimeGeneratedFunction(Meta.parse(_callbacks))
    cbset, ifelse_t0 = get_callbacks("https://xkcd.com/2694/") # Argument needed by @RuntimeGeneratedFunction
    compute_tstops = @RuntimeGeneratedFunction(Meta.parse(_tstops))

    # if model is written to file write the callback
    if write_to_file == true
        path_save = joinpath(dir_save, model_SBML.name * "_callbacks.jl")
        io = open(path_save, "w")
        write(io, _callbacks * "\n\n")
        write(io, _tstops)
        close(io)
    end

    if ret_all == true
        return reaction_system, specie_map, parameter_map, cbset, compute_tstops, ifelse_t0
    end

    # Depending on model return what is needed to perform forward simulations
    if isempty(model_SBML.events) && isempty(model_SBML.ifelse_bool_expressions)
        return reaction_system, specie_map, parameter_map
    end
    if isempty(model_SBML.ifelse_bool_expressions) && !isempty(model_SBML.events)
        verbose && @info "SBML model with events - output returned as odesys, specie_map, parameter_map, cbset, compute_tstops\nFor how to simulate model see documentation"
        return reaction_system, specie_map, parameter_map, cbset, compute_tstops
    end
    if isempty(model_SBML.ifelse_parameters) && !isempty(cbset)
        verbose && @info "SBML model with piecewise rewritten to event - output returned as odesys, specie_map, parameter_map, cbset, compute_tstops, cb_t0\nFor how to simulate model see documentation"
        return reaction_system, specie_map, parameter_map, cbset, compute_tstops, ifelse_t0
    end
end


function reactionsystem_from_SBML(model_SBML::ModelSBML)::String

    # Setup Catalyst ReactionNetwork
    _species_write, _species_write_array, _specie_map_write = SBMLImporter.get_specie_map(model_SBML, reaction_system=true)
    _parameters_write, _parameters_write_array, _parameter_map_write = SBMLImporter.get_parameter_map(model_SBML, reaction_system=true)

    # In case a specie (or parameter) appear as a rate-rule, algebraic or assignment rule they need to be treated as
    # MTK variable for the the downstream processing. This might turn the species block empty, then it must be removed
    _variables_write = "\tvs = ModelingToolkit.@variables"
    rule_variables = vcat(model_SBML.assignment_rule_variables, model_SBML.rate_rule_variables,
                          model_SBML.algebraic_rule_variables)
    filter!(x -> x ∉ keys(model_SBML.generated_ids), rule_variables)
    for variable in rule_variables
        _species_write = replace(_species_write, variable * "(t) " => "")
        _variables_write *= " " * variable * "(t)"
    end
    if !isempty(rule_variables)
        _species_write *= "\n" * _variables_write * "\n"
        no_species = all([any(occursin.(x, rule_variables)) for x in keys(model_SBML.species)])
    else
        no_species = false
    end

    # Reaction stoichiometry and propensities
    integer_stoichiometries::Bool = true
    _reactions_write = "\treactions = [\n"
    for (id, r) in model_SBML.reactions
        reactants_stoichiometries, reactants, integer_stoichiometries1 = get_reaction_side(r, :Reactants)
        products_stoichiometries, products, integer_stoichiometries2 = get_reaction_side(r, :Products)
        propensity = r.kinetic_math
        _reactions_write *= ("\t\tCatalyst.Reaction(" * propensity * ", " *
                             reactants * ", " * products * ", " *
                             reactants_stoichiometries * ", " * products_stoichiometries * "; only_use_rate=true),\n")

        # If it has already been assigned false we know that all stoichiometries are not
        # integer numbers, and if either of integer_stoichiometries are false integer_stoichiometries
        # should become false
        if integer_stoichiometries == true
            integer_stoichiometries = !any([integer_stoichiometries1, integer_stoichiometries2] .== false)
        end
    end

    # Rules are directly encoded into the Catalyst.Reaction vector
    for variable in vcat(model_SBML.rate_rule_variables, model_SBML.assignment_rule_variables)
        if variable ∈ keys(model_SBML.species)
            @unpack formula, assignment_rule = model_SBML.species[variable]
        elseif variable ∈ keys(model_SBML.parameters)
            @unpack formula, assignment_rule = model_SBML.parameters[variable]
        elseif variable ∈ keys(model_SBML.compartments)
            @unpack formula, assignment_rule = model_SBML.compartments[variable]
        else
            continue
        end
        if assignment_rule == true
            _reactions_write *= "\t\t" * variable * " ~ " * formula * ",\n"
        else
            _reactions_write *= "\t\tD(" * variable * ") ~ " * formula * ",\n"
        end
    end
    for formula in values(model_SBML.algebraic_rules)
        _reactions_write *= "\t\t" * formula * ",\n"
    end
    _reactions_write *= "\t]\n"

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
    _rn_write = "\trn = Catalyst.ReactionSystem(reactions, t, $sps_arg, ps; name=:" * model_SBML.name * ", combinatoric_ratelaws=$combinatoric_ratelaws_arg)"

    # Create a function returning the ReactionSystem, specie-map, and parameter-map
    _function_write = "function get_reaction_system(foo)\n\n"
    _function_write *= _species_write * "\n"
    _function_write *= _parameters_write * "\n\n"
    if !isempty(model_SBML.rate_rule_variables)
        _function_write *= "\tD = Differential(t)\n\n"
    end
    _function_write *= _reactions_write * "\n\n"
    _function_write *= _rn_write * "\n\n"
    _function_write *= _specie_map_write * "\n"
    _function_write *= _parameter_map_write * "\n"
    _function_write *= "\treturn rn, specie_map, parameter_map\nend"

    return _function_write
end



function get_reaction_side(r::ReactionSBML, which_side::Symbol)::Tuple{String, String, Bool}

    if which_side === :Reactants
        species, stoichiometries = r.reactants, r.reactants_stoichiometry
    elseif which_side === :Products
        species, stoichiometries = r.products, r.products_stoichiometry
    end

    # Case where we go from ϕ -> prod (or reverse)
    if isempty(species)
        return "nothing", "nothing", true
    end

    integer_stoichiometry::Bool = true
    _stoichiometries_str = "["
    for i in eachindex(species)
        stoichiometry, _integer_stoichiometry = parse_stoichiometry_reaction_system(stoichiometries[i])
        _stoichiometries_str *= stoichiometry * ", "
        if integer_stoichiometry == true
            integer_stoichiometry = _integer_stoichiometry
        end
    end
    _stoichiometries_str = _stoichiometries_str[1:end-2] * "]"
    _species_str = "[" * prod([s * ", " for s in species])[1:end-2] * "]"

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