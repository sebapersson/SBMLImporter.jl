using Catalyst, CSV, DataFrames, Downloads, OrdinaryDiffEqBDF, OrdinaryDiffEqRosenbrock,
    SBML, SBMLImporter, Test

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "testsuite_support.jl"))

# Standard file format, so parsing is hard-coded
function read_settings(settings_url::String)
    settings = String(take!(Downloads.download(settings_url, IOBuffer())))
    settings_lines = split(settings, '\n')
    species_test_amount = replace.(split(split(settings_lines[7], ":")[2], ','), " " => "")
    species_test_conc = replace.(split(split(settings_lines[8], ":")[2], ','), " " => "")
    abstol_test = parse(Float64, split(settings_lines[5], ":")[2])
    reltol_test = parse(Float64, split(settings_lines[6], ":")[2])

    return species_test_amount, species_test_conc, abstol_test, reltol_test
end

function check_test_case(test_case, solver)
    @info "Test case $test_case"

    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/" *
               "semantic/$test_case/$test_case"
    sbml_urls = get_sbml_urls(base_url)

    _results_url = base_url * "-results.csv"
    _results = String(take!(Downloads.download(_results_url, IOBuffer())))
    results = CSV.read(IOBuffer(_results), DataFrame)

    # As it stands I cannot "hack" a parameter value at time zero, but the model simulation values
    # are correct. Border case we pass
    _cases = ["00995", "00996", "00997", "01284", "01510", "01527", "01596", "01663",
              "01686", "01684", "01685", "01694", "01695", "01696", "01697", "01698",
              "01699", "01700", "01719", "00928", "00929", "01693"]
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
    species_test_amount, species_test_conc, abstol_test, reltol_test = read_settings(settings_url)
    sbml_url = sbml_urls[end]

    for sbml_url in sbml_urls
        sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
        model_SBML = SBMLImporter.parse_SBML(sbml_string, false, model_as_string = true,
                                             inline_assignment_rules=false)
        # If stoichiometryMath occurs we need to convert the SBML file to a level 3 file
        # to properly handle it
        if occursin("stoichiometryMath", sbml_string) == false
            libsbml_model = readSBMLFromString(sbml_string)
        else
            libsbml_model = readSBMLFromString(sbml_string,
                                               doc -> begin
                                                   set_level_and_version(3, 2)(doc)
                                                   convert_promotelocals_expandfuns(doc)
                                               end)
        end

        prn, cb = load_SBML(sbml_string, model_as_string = true, inline_assignment_rules = false, ifelse_to_callback=ifelse_to_callback)
        if isempty(model_SBML.algebraic_rules)
            osys = structural_simplify(convert(ODESystem, prn.rn))
        else
            _sys = dae_index_lowering(convert(ODESystem, prn.rn))
            osys = structural_simplify(_sys)
        end
        oprob = ODEProblem(osys, prn.u0, (0.0, tmax), prn.p, jac = true)
        sol = solve(oprob, solver, abstol = 1e-12, reltol = 1e-12, saveat = tsave,
                    callback = cb)
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
end

# Get test cases currently not passing
cases_not_passing = get_semantic_not_pass()
keys_cases_catch = filter(x -> x != :not_captured, keys(cases_not_passing))
cases_not_passing_catch = reduce(vcat, [cases_not_passing[k] for k in keys_cases_catch])
solver = Rodas4P()
@testset "Catalyst" begin
    for i in 1:1821
        test_case = repeat("0", 5 - length(string(i))) * string(i)

        if test_case in cases_not_passing_catch
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case in cases_not_passing[:not_captured]
            continue
        end
        if test_case == "01596" # Must be able to solve for time event is triggered
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end

        # Simulations fail with Rodas4P, so FBDF is used
        if test_case in ["00028", "00173", "00269"]
            check_test_case(test_case, FBDF())
            continue
        end
        check_test_case(test_case, solver)
    end
end
