function _get_odeproblem(
        model_SBML_prob::ModelSBMLProb, model_SBML::ModelSBML
    )::Tuple{SciMLBase.ODEProblem, Function, CallbackSet}
    _fode = _template_odeproblem(model_SBML_prob, model_SBML)
    fode = @RuntimeGeneratedFunction(Meta.parse(_fode))

    # Parameter vector is returned as a ComponentArray
    psvals = parse.(Float64, last.(model_SBML_prob.psmap))
    psvec = ComponentArray(; (Symbol.(first.(model_SBML_prob.psmap)) .=> psvals)...)

    # Initial values returned as a function u0(p, t0) to correctly map initial values
    _fu0 = _template_u0_odeproblem(model_SBML_prob, model_SBML)
    fu0 = @RuntimeGeneratedFunction(Meta.parse(_fu0))
    u0tmp = zeros(Float64, length(model_SBML_prob.umodel))
    oprob = SciMLBase.ODEProblem(fode, u0tmp, (0.0, 10.0), psvec)

    # Potential model callbacks. p_PEtab and _specie_ids are needed to inform the importer
    # about parameter and specie ordering in the ODEProblem. This information is not
    # needed if the model is parsed into a System
    cb = create_callbacks(
        oprob, model_SBML, model_SBML.name; p_PEtab = model_SBML_prob.ps,
        float_tspan = true, _specie_ids = model_SBML_prob.umodel
    )
    return oprob, fu0, cb
end

function ModelSBMLProb(model_SBML)
    species, variables, specie_map = _get_system_variables(model_SBML, true)
    umodel = _get_umodel(species, variables)
    umap = _get_ps_or_umap(specie_map, :u)
    @assert length(umap) == length(umodel) "Map and u length does not match in ODE parsing"

    parameters, parameter_map = _get_system_parameters(model_SBML)
    ps = split(parameters[28:(end - 1)], " ") .|> string
    psmap = _get_ps_or_umap(parameter_map, :ps)
    @assert length(ps) == length(psmap) "Map and ps length does not match in ODE parsing"

    # RHS ODEs (in the order of umodel)
    odes = fill("", length(umodel))
    for (i, uid) in pairs(umodel)
        if isempty(model_SBML.species[uid].formula)
            formula = "0.0"
        else
            formula = model_SBML.species[uid].formula
        end
        odes[i] = "du[$i] = " * formula * "\n"
    end

    return ModelSBMLProb(umodel, umap, ps, psmap, odes)
end

function _get_umodel(species::String, variables::String)::Vector{String}
    umodel = split(species[26:(end - 1)] * variables[33:(end - 1)], " ") .|> string
    return replace.(umodel, "(t)" => "")
end

function _get_ps_or_umap(map::String, which::Symbol)::Vector{Pair{String, String}}
    if which == :u
        _map = split(map[16:(end - 1)], ",\n\t")[1:(end - 1)]
    else
        _map = split(map[19:(end - 1)], ",\n\t")[1:(end - 1)]
    end
    _map = replace.(_map, "\n\t" => "")
    _map = split.(_map, " =>")
    _map = (first.(_map) .=> last.(_map))
    return _map
end
