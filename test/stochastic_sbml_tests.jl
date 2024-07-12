using SBMLImporter
using CSV
using DataFrames
using OrdinaryDiffEq
using ModelingToolkit
using SBML
using Test
using Downloads
using Catalyst
using JumpProcesses

include(joinpath(@__DIR__, "common.jl"))

function read_settings(settings_url::String)
    settings = String(take!(Downloads.download(settings_url, IOBuffer())))
    settings_lines = split(settings, '\n')
    species_test = Symbol.(replace.(split(split(settings_lines[4], ":")[2], ','),
                                    " " => ""))

    mean_range = parse.(Float64,
                        split(replace.(split(settings_lines[10], ":")[2], r"\(|\)| " => ""),
                              ","))
    sd_range = parse.(Float64,
                      split(replace.(split(settings_lines[11], ":")[2], r"\(|\)| " => ""),
                            ","))

    return species_test, mean_range, sd_range
end

function test_stochastic_testcase(test_case::String; nsolve::Integer = 20000)
    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/stochastic/$test_case/$test_case"

    species_test, mean_range, sd_range = read_settings(base_url * "-settings.txt")

    _results_url = base_url * "-results.csv"
    _results = String(take!(Downloads.download(_results_url, IOBuffer())))
    results = CSV.read(IOBuffer(_results), DataFrame)
    t_save = "Time" in names(results) ? Float64.(results[!, :Time]) :
             Float64.(results[!, :time])
    t_save = Vector{Float64}(t_save)
    tmax = maximum(t_save)

    sbml_urls = get_sbml_urls(base_url)
    sbml_url = sbml_urls[end]

    # To reduce runtime for the most demanding test-case
    if test_case == "00005"
        sbml_urls = [sbml_urls[end]]
    end

    for sbml_url in sbml_urls
        sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
        # c = n / V => n = c * V

        model_SBML = SBMLImporter.build_SBML_model(sbml_string, model_as_string = true)
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

        #SBMLImporter.SBML_to_ODESystem(sbml_string, ret_all=true, model_as_string=true, write_to_file=true)
        parsed_rn, cb = load_SBML(sbml_string, model_as_string = true,
                                  inline_assignment_rules = false)
        tspan = (0.0, tmax)
        dprob = DiscreteProblem(parsed_rn.rn, parsed_rn.u₀, tspan, parsed_rn.p)
        jprob = JumpProblem(parsed_rn.rn, dprob, Direct(); save_positions = (false, false))
        eprob = EnsembleProblem(jprob)
        if test_case != "00033"
            sol = solve(eprob, SSAStepper(), EnsembleSerial(), trajectories = nsolve,
                        saveat = t_save, callback = cb)
        else
            sol = solve(eprob, FunctionMap(), EnsembleSerial(), trajectories = nsolve,
                        saveat = t_save, callback = cb)
        end

        for to_check in species_test
            to_check_no_whitespace = replace(string(to_check), " " => "")
            i_specie = findfirst(x -> x == to_check_no_whitespace,
                                 replace.(string.(states(jprob.prob.f.sys)), "(t)" => ""))
            sim_mean, sim_var = SciMLBase.EnsembleAnalysis.timepoint_meanvar(sol, t_save)
            sim_mean = sim_mean[i_specie, :]
            sim_sd = sqrt.(sim_var[i_specie, :])

            if to_check_no_whitespace * "-mean" ∈ names(results)
                reference_sol = results[!, to_check_no_whitespace * "-mean"]
                @test all(reference_sol .+ mean_range[1] .< sim_mean .<
                          reference_sol .+ mean_range[2])
            end

            if to_check_no_whitespace * "-sd" ∈ names(results)
                reference_sol = results[!, to_check_no_whitespace * "-sd"]
                @test all(reference_sol .+ sd_range[1] .< sim_sd .<
                          reference_sol .+ sd_range[2])
            end
        end
    end
end

@testset "Stochastic tests" begin
    for i in 1:39

        # 19 has assignment rule
        # 33 has continous callback
        if i ∈ [19, 33]
            continue
        end

        @info "Test case $i"
        test_case = repeat("0", 5 - length(string(i))) * string(i)

        # i = 5 nast test case were many simulations must be run to get a
        # good mean estimate (as amounts are large)
        if i == 5
            test_stochastic_testcase(test_case; nsolve = 70000)
        else
            test_stochastic_testcase(test_case)
        end
    end
end
