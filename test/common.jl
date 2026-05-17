function get_sbml_urls(base_url::String)
    levels = ["2", "3"]
    sublevels = ["1", "2", "3", "4", "5"]
    sbml_urls = String[]
    for level in levels
        for sublevel in sublevels
            sbml_url = base_url * "-sbml-l" * level * "v" * sublevel * ".xml"
            try
                sbml = String(take!(Downloads.download(sbml_url, IOBuffer())))
                push!(sbml_urls, sbml_url)
            catch
            end
        end
    end
    return sbml_urls
end

# Function to get model-str
function get_model_str(test_case)
    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/semantic/$test_case/$test_case"
    sbml_urls = get_sbml_urls(base_url)
    sbml_url = sbml_urls[end]
    sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
    model_SBML = SBMLImporter.build_SBML_model(
        sbml_string, model_as_string = true, inline_assignment_rules = false
    )

    parsed_model_SBML = SBMLImporter._reactionsystem_from_SBML(model_SBML)
    model_str = SBMLImporter.reactionsystem_to_string(
        parsed_model_SBML, false, "", model_SBML
    )

    return model_str, model_SBML
end

function get_test_tags(test_case)
    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/semantic/$test_case/$test_case"
    test_file_url = base_url * "-model.m"
    test_file_str = String(take!(Downloads.download(test_file_url, IOBuffer())))
    return extract_tags(test_file_str)
end

function extract_tags(txt::AbstractString)
    pat = r"^\s*(componentTags|testTags)\s*:\s*(.*)$"m
    tags = Dict{Symbol, Vector{String}}()
    for m in eachmatch(pat, txt)
        key = m.captures[1]
        value = strip.(split(m.captures[2], ','))
        tags[Symbol(key)] = value
    end
    return tags
end

# Standard file format, so parsing is hard-coded
function read_semantic_settings(settings_url::String)
    settings = String(take!(Downloads.download(settings_url, IOBuffer())))
    settings_lines = split(settings, '\n')
    species_test_amount = replace.(split(split(settings_lines[7], ":")[2], ','), " " => "")
    species_test_conc = replace.(split(split(settings_lines[8], ":")[2], ','), " " => "")
    abstol_test = parse(Float64, split(settings_lines[5], ":")[2])
    reltol_test = parse(Float64, split(settings_lines[6], ":")[2])

    return species_test_amount, species_test_conc, abstol_test, reltol_test
end

function check_semantic_test_case(test_case, solver)
    @info "Test case $test_case"

    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/" *
        "semantic/$test_case/$test_case"
    sbml_urls = get_sbml_urls(base_url)

    _results_url = base_url * "-results.csv"
    _results = String(take!(Downloads.download(_results_url, IOBuffer())))
    results = CSV.read(IOBuffer(_results), DataFrame)

    # As it stands I cannot "hack" a parameter value at time zero, but the model simulation values
    # are correct. Border case we pass
    _cases = [
        "00995", "00996", "00997", "01284", "01510", "01527", "01596", "01663",
        "01686", "01684", "01685", "01694", "01695", "01696", "01697", "01698",
        "01699", "01700", "01719", "00928", "00929", "01693",
    ]
    if test_case in _cases
        results = results[2:end, :]
    end
    # Have trigger time t = 0, which is not currently handled
    if test_case == "01488"
        ifelse_to_callback = false
    else
        ifelse_to_callback = true
    end
    colnames = Symbol.(replace.(string.(names(results)), " " => ""))
    rename!(results, colnames)

    # Case for FBA models an exceptions should be thrown
    if !("Time" in names(results) || "time" in names(results))
        sbml_string = String(take!(Downloads.download(sbml_urls[1], IOBuffer())))
        SBMLImporter.parse_SBML(sbml_string, false; model_as_string = true)
    end

    tsave = "Time" in names(results) ? Float64.(results[!, :Time]) : Float64.(results[!, :time])
    tsave = Vector{Float64}(tsave)
    tmax = maximum(tsave)
    components_test = filter(x -> x ∉ [:time, :Time], colnames)

    # Read settings file
    settings_url = base_url * "-settings.txt"
    species_test_amount, species_test_conc, abstol_test, reltol_test = read_semantic_settings(settings_url)

    for sbml_url in sbml_urls
        sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
        model_SBML = SBMLImporter.parse_SBML(
            sbml_string, false, model_as_string = true, inline_assignment_rules = false
        )
        # If stoichiometryMath occurs we need to convert the SBML file to a level 3 file
        # to properly handle it
        if occursin("stoichiometryMath", sbml_string) == false
            libsbml_model = readSBMLFromString(sbml_string)
        else
            libsbml_model = readSBMLFromString(
                sbml_string,
                doc -> begin
                    set_level_and_version(3, 2)(doc)
                    convert_promotelocals_expandfuns(doc)
                end
            )
        end

        rn, cb = load_SBML(
            sbml_string, model_as_string = true, inline_assignment_rules = false,
            ifelse_to_callback = ifelse_to_callback
        )
        u0 = get_u0_map(rn)
        ps = get_parameter_map(rn)
        if isempty(model_SBML.algebraic_rules)
            osys = ModelingToolkitBase.mtkcompile(ode_model(rn))
        else
            _sys = ModelingToolkit.dae_index_lowering(ode_model(rn))
            osys = ModelingToolkitBase.mtkcompile(_sys)
        end
        oprob = ODEProblem(osys, merge(Dict(u0), Dict(ps)), (0.0, tmax), jac = true)
        sol = solve(
            oprob, solver, abstol = 1.0e-12, reltol = 1.0e-12, saveat = tsave,
            callback = cb,
        )
        model_parameters = string.(parameters(sol.prob.f.sys))

        for _component_test in components_test
            component_test = replace(string(_component_test), " " => "")
            expected_res = results[!, Symbol(component_test)]

            if component_test in model_parameters
                p = sol.prob.ps[Symbol(component_test)]
                result_is_inf = all(isinf.(expected_res))
                result_is_nan = all(isnan.(expected_res))
                result_is_pos = all(expected_res .> 0)

                if result_is_inf && result_is_pos
                    @test isinf(p) && sol.prob.p[ip] > 0
                elseif result_is_inf && !result_is_pos
                    @test isinf(p) && sol.prob.p[ip] < 0
                elseif result_is_nan
                    @test isnan(p)
                else
                    absdiff = abs.(p .- expected_res)
                    @test all(absdiff .< abstol_test .+ reltol_test .* abs.(expected_res))
                end
                continue
            end

            is_sbml_specie = haskey(libsbml_model.species, component_test)
            if component_test in vcat(species_test_conc, species_test_amount) && is_sbml_specie
                compartment_name = libsbml_model.species[component_test].compartment
                unit = model_SBML.species[component_test].unit
                if unit == :Concentration && component_test ∉ species_test_amount
                    c = 1.0
                elseif compartment_name in model_parameters
                    c = sol.prob.ps[Symbol(compartment_name)]
                else
                    c = sol[Symbol(compartment_name)]
                end
                if component_test in species_test_amount && unit == :Amount
                    c = 1.0
                elseif component_test in species_test_amount
                    c = 1 ./ c
                end
            else
                c = 1.0
            end
            absdiff = abs.(sol[Symbol(component_test)] ./ c .- expected_res)
            @test all(absdiff .< abstol_test .+ reltol_test .* abs.(expected_res))
        end
    end
    return nothing
end
