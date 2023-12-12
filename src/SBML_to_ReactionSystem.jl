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

    # Things I am pretty certain we do not support yet
    if !isempty(model_SBML.algebraic_rules)
        throw(SBMLSupport("For ReactionSystem we do not support algebraic rules"))
    end
    if !isempty(model_SBML.rate_rule_variables)
        throw(SBMLSupport("For ReactionSystem we do not support rate rules"))
    end
    # Assignment rules might be handled correctly, however, I am not willing 
    # to die on the hill, but will be fixed downstream 

    # Build the ReactionSystem
    _model = reactionsystem_from_SBML(model_SBML)
    _get_reaction_system = @RuntimeGeneratedFunction(Meta.parse(_model))
    reaction_system, specie_map, parameter_map = _get_reaction_system("https://xkcd.com/303/") # Argument needed by @RuntimeGeneratedFunction

    # Build callback functions
    _callbacks, _tstops = create_callbacks_SBML(ode_system, model_SBML, model_SBML.name)
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

    # Reaction stoichiometry and propensities
    _reactions_write = "\treactions = [\n"
    for (id, r) in model_SBML.reactions
        reactants_stoichiometries, reactants = get_reaction_side(r, :Reactants)
        products_stoichiometries, products = get_reaction_side(r, :Products)
        propensity = r.kinetic_math
        _reactions_write *= ("\t\tCatalyst.Reaction(" * propensity * ", " * 
                             reactants * ", " * products * ", " *
                             reactants_stoichiometries * ", " * products_stoichiometries * "; only_use_rate=true),\n")
    end
    _reactions_write *= "\t]\n"
    
    # ReactionSystem
    _rn_write = "\trn = Catalyst.ReactionSystem(reactions, t, sps, ps; name=:" * model_SBML.name * ")"

    # Create a function returning the ReactionSystem, specie-map, and parameter-map
    _function_write = "function get_reaction_system(foo)\n\n"
    _function_write *= _species_write * "\n"
    _function_write *= _parameters_write * "\n\n"
    _function_write *= _reactions_write * "\n\n"
    _function_write *= _rn_write * "\n\n"
    _function_write *= _specie_map_write * "\n"
    _function_write *= _parameter_map_write * "\n"
    _function_write *= "\treturn rn, specie_map, parameter_map\nend"

    return _function_write
end



function get_reaction_side(r::ReactionSBML, which_side::Symbol)::Tuple{String, String}

    if which_side === :Reactants
        species, stoichiometries = r.reactants, r.reactants_stoichiometry
    elseif which_side === :Products
        species, stoichiometries = r.products, r.products_stoichiometry
    end

    _stoichiometries_str = "["
    for i in eachindex(species)
        stoichiometry = parse_stoichiometry_reaction_system(stoichiometries[i])
        _stoichiometries_str *= stoichiometry * ", "
    end
    _stoichiometries_str = _stoichiometries_str[1:end-2] * "]"
    _species_str = "[" * prod([s * ", " for s in species])[1:end-2] * "]"

    if isempty(_stoichiometries_str)
        return "nothing", "nothing"
    else
        return _stoichiometries_str, _species_str
    end
end


# If possible parse stoichiometry to an integer
function parse_stoichiometry_reaction_system(stoichiometry::String)::String
    if isnothing(tryparse(Float64, stoichiometry))
        return stoichiometry[i] * "*"
    end
    _stoichiometry = tryparse(Float64, stoichiometry)
    try
        return string(Int64(_stoichiometry))
    catch
        return string(_stoichiometry)
    end
end