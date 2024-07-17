function _get_reaction_system(model_SBML_sys::ModelSBMLSystem, name::String)
    # The ReactionSystem must be built via eval, as creating a function that returns
    # the rn fails for large models
    eval(Meta.parse("ModelingToolkit.@variables t"))
    eval(Meta.parse("D = Catalyst.Differential(t)"))
    # A model must have either variables or species, which dictates the sps call to
    # the reaction system
    if model_SBML_sys.has_species == true
        sps = eval(Meta.parse(model_SBML_sys.species))
    else
        sps = Any[]
    end
    if model_SBML_sys.variables != "\tvs = ModelingToolkit.@variables "
        vs = eval(Meta.parse(model_SBML_sys.variables))
    else
        vs = Any[]
    end
    if vs == Any[]
        sps_arg = sps
    elseif model_SBML_sys.has_species == true
        sps_arg = [sps; vs]
    else
        sps_arg = vs
    end

    if model_SBML_sys.parameters != "\tps = Catalyst.@parameters "
        ps = eval(Meta.parse(model_SBML_sys.parameters))
    else
        ps = Any[]
    end

    # Build reaction system from its components. The reactions must be parsed sequentially
    # to avoid any stackoverflow error for very large models with more than 10^4 reactions
    # TODO: This is going to be a pain with with Catalyst v14... (maybe can then add vector)
    reactions = _reaction_str_to_vec(model_SBML_sys)
    r1 = eval(Meta.parse(reactions[1]))
    r1 = isnothing(r1) ? Catalyst.Reaction[] : [r1]
    rn = Catalyst.ReactionSystem(r1, t, sps_arg, ps; name = Symbol(name),
                                 combinatoric_ratelaws = model_SBML_sys.all_integer_S)
    for i in eachindex(reactions)
        i == 1 && continue
        # Adding a reaction and equation differs, as with reaction we can simply add via
        # addreaction! (much easier)
        r = eval(Meta.parse(reactions[i]))
        if occursin("Reaction", reactions[i])
            Catalyst.addreaction!(rn, r)
        else
            Catalyst.reset_networkproperties!(rn)
            push!(ModelingToolkit.get_eqs(rn), r)
            sort(ModelingToolkit.get_eqs(rn); by = Catalyst.eqsortby)
        end
    end
    specie_map = eval(Meta.parse(model_SBML_sys.specie_map))
    parameter_map = eval(Meta.parse(model_SBML_sys.parameter_map))
    return rn, specie_map, parameter_map
end

function write_reactionsystem(model_SBML_sys::ModelSBMLSystem, dirsave::String,
                              model_SBML::ModelSBML)::String
    # If model is written to file save it in the same directory as the SBML-file. Only
    # save if model is not provided as a string (as then there is not file)
    pathsave = joinpath(dirsave, model_SBML.name * ".jl")

    sps = model_SBML_sys.has_species ? model_SBML_sys.species : "Any[]"
    if model_SBML_sys.variables != "\tvs = ModelingToolkit.@variables "
        vs = model_SBML_sys.variables
    else
        vs = "Any[]"
    end
    if vs == "Any[]"
        sps_arg = "sps"
    elseif model_SBML_sys.has_species == true
        sps_arg = "[sps; vs]"
    else
        sps_arg = "vs"
    end
    if model_SBML_sys.parameters != "\tps = Catalyst.@parameters "
        ps = model_SBML_sys.parameters
    else
        ps = "Any[]"
    end
    reactions = model_SBML_sys.reactions
    model_name = "\"" * model_SBML.name * "\""
    comb_ratelaws = string(model_SBML_sys.all_integer_S)
    rn_write = "Catalyst.ReactionSystem(_reactions, t, sps_arg, ps; " *
               "name=Symbol($(model_name)), " *
               "combinatoric_ratelaws=comb_ratelaws)"

    frn = "function get_reaction_system(foo)\n"
    frn *= "\tModelingToolkit.@variables t\n"
    frn *= "\tD = Differential(t)\n"
    frn *= sps * "\n"
    frn *= vs * "\n"
    frn *= "\tsps_arg = " * sps_arg * "\n"
    frn *= ps * "\n\n"
    frn *= reactions * "\n"
    frn *= "\tcomb_ratelaws = " * comb_ratelaws * "\n"
    frn *= "\trn = " * rn_write * "\n\n"
    frn *= model_SBML_sys.specie_map * "\n\n"
    frn *= model_SBML_sys.parameter_map * "\n"
    frn *= "\treturn rn, specie_map, parameter_map\nend"

    open(pathsave, "w") do f
        write(f, frn)
    end
    return frn
end

function _to_system_syntax(model_SBML::ModelSBML, inline_assignment_rules::Bool,
                           massaction::Bool)::ModelSBMLSystem
    # If model is empty of derivatives a dummy state must be added to be able to create
    # a ReactionSystem
    if _has_derivatives(model_SBML) == false
        model_SBML.species["__foo__"] = SpecieSBML("__foo__", false, false, "1.0", "0.0",
                                                   "1.0", "", :Amount, false, false, false,
                                                   false, false, false, false)
    end

    species, variables, specie_map = _get_system_variables(model_SBML,
                                                           inline_assignment_rules)
    parameters, parameter_map = _get_system_parameters(model_SBML)

    reactions, all_integer_S = _get_system_reactions(model_SBML, massaction)

    # Rules (ODEs) go into the Catalyst reaction vector when parsing the system. Note, we
    # need to close with "\t]\n" for nice printing to file
    reactions *= _get_system_rules(model_SBML)
    reactions *= "\t]\n"

    has_species = species != "\tsps = Catalyst.@species "
    return ModelSBMLSystem(species, specie_map, variables, parameters, parameter_map,
                           reactions, has_species, all_integer_S)
end

function _get_system_variables(model_SBML::ModelSBML,
                               inline_assignment_rules::Bool)::Tuple{String, String, String}
    # Species in reactions are treated as Catalyst.@species, meanwhile, dynamic variables
    # not changed by reactions (rule variables) are ModelingToolkit variables. Note,
    # both mtk_variables and catalyst_species initial values are set via the specie_map
    mtk_variables = "\tvs = ModelingToolkit.@variables "
    catalyst_species = "\tsps = Catalyst.@species "
    specie_map = "\tspecie_map = [\n"
    for (specie_id, specie) in model_SBML.species
        specie_id in model_SBML.rule_variables && continue
        catalyst_species *= specie_id * "(t) "
        specie_map *= _template_value_map(specie_id, specie.initial_value)
    end
    for variable_id in model_SBML.rule_variables
        variable = _get_model_variable(variable_id, model_SBML)
        # If assignment rule is inlined, it is no longer a model variable, and thus no
        # longer treated as a model variable
        if variable.assignment_rule && inline_assignment_rules
            continue
        end
        if !(variable isa SpecieSBML) && variable.constant == true && continue
            continue
        end
        mtk_variables *= variable_id * "(t) "
        if !(variable_id in model_SBML.algebraic_rule_variables)
            specie_map *= _template_value_map(variable_id, variable.initial_value)
        end
    end
    specie_map *= "\t]"
    return catalyst_species, mtk_variables, specie_map
end

function _get_system_parameters(model_SBML::ModelSBML)::Tuple{String, String}
    catalyst_parameters = "\tps = Catalyst.@parameters "
    parameter_map = "\tparameter_map = [\n"
    variables = Iterators.flatten((model_SBML.parameters, model_SBML.compartments))
    for (variable_id, variable) in variables
        variable.constant == false && continue
        variable.assignment_rule == true && continue
        catalyst_parameters *= variable_id * " "
        parameter_map *= _template_value_map(variable_id, variable.formula)
    end
    parameter_map *= "\t]"
    return catalyst_parameters, parameter_map
end

function _get_system_reactions(model_SBML::ModelSBML, massaction::Bool)::Tuple{String, Bool}
    all_integer_S::Bool = true
    catalyst_reactions = "\t_reactions = [\n"
    for (id, r) in model_SBML.reactions
        # Can happen for models with species that are boundary conditions, as boundary
        # condition species do not partake in reactions
        if all(r.products .== "nothing") && all(r.reactants .== "nothing")
            continue
        end

        r_S, reactants, r_integer_S = _get_reaction_system_info(r, model_SBML, :Reactants)
        p_S, products, p_integer_S = _get_reaction_system_info(r, model_SBML, :Products)
        propensity = r.kinetic_math

        # A catalyst reaction is massaction (and can be efficiently simulated) if:
        # 1. Involved species has only_substance_unit=true
        # 2. Propensity does not include time t, or depend on rate-rule or assignment rules
        # 3. Stoichiometry is mass-action, e.g. stoichiometry is not given by model
        #    parameters
        is_massaction = _is_massaction(r, model_SBML, r_integer_S, p_integer_S)
        if massaction == true && is_massaction == false
            @warn "That the system is massaction was provided, however, the reaction " *
                  "$id is not massaction. It is parsed as non massaction reaction, which " *
                  "can negative impact simulation times."
        end
        catalyst_reactions *= _template_reaction(reactants, products, r_S, p_S, propensity,
                                                 is_massaction)

        # For creating the reaction system we need to know if only integer stoichiometry
        # holds for all reactions
        all_integer_S = all_integer_S ? all([r_integer_S, p_integer_S]) : all_integer_S
    end
    return catalyst_reactions, all_integer_S
end

function _get_system_rules(model_SBML::ModelSBML)::String
    catalyst_rules = ""
    rule_ids = unique(Iterators.flatten((model_SBML.rate_rule_variables,
                                         model_SBML.assignment_rule_variables)))
    for rule_id in rule_ids
        variable = _get_model_variable(rule_id, model_SBML)
        if rule_id in model_SBML.rate_rule_variables
            catalyst_rules *= _template_rate_rule(variable.name, variable.formula)
        else
            catalyst_rules *= _template_assignment_rule(variable.name, variable.formula)
        end
    end
    # Algebraic rules do not (as above) have an associated id
    for rule in values(model_SBML.algebraic_rules)
        catalyst_rules *= "\t\t" * rule.formula * ",\n"
    end
    return catalyst_rules
end

function _has_derivatives(model_SBML::ModelSBML)::Bool
    @unpack species, parameters, compartments = model_SBML
    any_non_assignment_s = any([s.assignment_rule == false for s in values(species)])
    any_rate_p = any([p.rate_rule for p in values(parameters)])
    any_rate_c = any([c.rate_rule for c in values(compartments)])
    if any_non_assignment_s == false && any_rate_c == false && any_rate_p == false
        return false
    else
        return true
    end
end

function _get_reaction_stoichiometry(stoichiometry::String)::Tuple{String, Bool}
    # stoichiometry can be a parameter or float, but, if it is an Int this means that
    # we can treat the reaction as an efficient mass-action reaction
    if !isnothing(tryparse(Float64, stoichiometry))
        _stoichiometry = parse(Float64, stoichiometry)
        try
            return string(Int64(_stoichiometry)), true
        catch
            return string(_stoichiometry), false
        end
    end
    return stoichiometry, false
end

function _get_reaction_system_info(r::ReactionSBML, model_SBML::ModelSBML,
                                   which_side::Symbol)::Tuple{String, String, Bool}
    if which_side === :Reactants
        species_id, stoichiometries = r.reactants, r.reactants_stoichiometry
    elseif which_side === :Products
        species_id, stoichiometries = r.products, r.products_stoichiometry
    end

    # When going from ϕ -> product, or reactant -> ϕ
    # Second condition can happen for boundary conditions
    if isempty(species_id) || all(species_id .== "nothing")
        return "nothing", "nothing", true
    end

    # In a SBML.Reaction double counting can occur (same specie appears twice), this is not
    # allowed by Catalyst, and must be considered when parsing stoichiometries (referred to
    # as S below).
    # If stoichiometries are integers, the system can be simulated with efficient massaction
    # simulators, hence, integer_S is checked
    integer_S::Bool = true
    S = fill("", length(filter(x -> x != "nothing", unique(species_id))))
    species_parsed = fill("", length(S))
    k = 1
    for (i, specie_id) in pairs(species_id)
        # Boundary condition, should not affect model dynamics
        if specie_id == "nothing"
            continue
        end

        # SBML models can have conversion factors that scale stoichiometry, if this
        # is the case, the stoichiometry is no longer an integer
        _S, _integer_S = _get_reaction_stoichiometry(stoichiometries[i])
        cf_scaling = _get_cf_scaling(model_SBML.species[specie_id], model_SBML)
        _S = isempty(cf_scaling) ? _S : _apply(*, _S, cf_scaling)
        _integer_S *= isempty(cf_scaling)

        if specie_id in species_parsed
            is = findfirst(x -> x == specie_id, species_parsed)
            S[is] = _apply(+, S[is], _S)
        else
            species_parsed[k] = specie_id
            S[k] = _S
            k += 1
        end
        integer_S = integer_S == true ? _integer_S : integer_S
    end
    catalyst_S = "[" * prod([_S * ", " for _S in S])[1:(end - 2)] * "]"
    catalyst_species = "[" * prod([s * ", " for s in species_parsed])[1:(end - 2)] * "]"
    return catalyst_S, catalyst_species, integer_S
end

function _is_massaction(r::ReactionSBML, model_SBML::ModelSBML, r_integer_S::Bool,
                        p_integer_S::Bool)::Bool
    r_integer_S == false && return false
    p_integer_S == false && return false
    r.stoichiometry_mass_action == false && return false
    for specie in Iterators.flatten((r.reactants, r.products))
        specie == "nothing" && continue
        if model_SBML.species[specie].only_substance_units == false
            return false
        end
    end
    formula = r.kinetic_math
    for rule_variable in model_SBML.rule_variables
        if _replace_variable(formula, rule_variable, "") != formula
            return false
        end
    end
    return true
end

function _reaction_str_to_vec(model_SBML_sys::ModelSBMLSystem)
    _reactions = model_SBML_sys.reactions[15:(end - 5)]
    _reactions = replace(_reactions, '\t' => "")
    _reactions = split(_reactions, ",\n")
    _reactions[1] = _reactions[1][3:end]
    return _reactions
end

function update_rate_reaction(rx; combinatoric_ratelaw::Bool = true)
    @set rx.rate = Catalyst.simplify((rx.rate^2) /
                                     Catalyst.oderatelaw(rx;
                                                         combinatoric_ratelaw = combinatoric_ratelaw))
end

function _get_dir_save(write_to_file::Bool, model_as_string::Bool,
                       path::String)::Union{Nothing, String}
    if write_to_file == true
        dir_save = if model_as_string
            nothing
        else
            joinpath(splitdir(path)[1], "SBML")
        end
    else
        dir_save = nothing
    end
    if !(isnothing(dir_save) || isdir(dir_save))
        mkdir(dir_save)
    end
    return dir_save
end
