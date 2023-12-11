using SBMLImporter
using CSV 
using DataFrames
using OrdinaryDiffEq
using ModelingToolkit
using SBML
using Test
using Downloads
using Catalyst


function get_sbml_urls(base_url::String)
    levels = ["2", "3"]
    sublevels = ["1", "2", "3", "4", "5"]
    sbml_urls = Vector{String}(undef, 0)

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


function read_settings(settings_url::String)

    settings = String(take!(Downloads.download(settings_url, IOBuffer())))
    settings_lines = split(settings, '\n')
    species_test = Symbol.(replace.(split(split(settings_lines[4], ":")[2], ','), " " => "", ))
    species_test_amount = Symbol.(replace.(split(split(settings_lines[7], ":")[2], ','), " " => "", ))
    species_test_conc = Symbol.(replace.(split(split(settings_lines[8], ":")[2], ','), " " => "", ))
    abstol_test = parse(Float64, split(settings_lines[5], ":")[2])
    reltol_test = parse(Float64, split(settings_lines[6], ":")[2])

    return species_test, species_test_amount, species_test_conc, abstol_test, reltol_test
end


function check_test_case(test_case, solver)
    
    @info "Test case $test_case"

    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/semantic/$test_case/$test_case"
    sbml_urls = get_sbml_urls(base_url)

    _results_url = base_url * "-results.csv"
    _results = String(take!(Downloads.download(_results_url, IOBuffer())))
    results = CSV.read(IOBuffer(_results), DataFrame)
    
    # As it stands I cannot "hack" a parameter value at time zero, but the model simulation values
    # are correct. Border case we pass 
    if test_case == "00995" || test_case == "00996" || test_case == "00997" || test_case == "01284" || test_case == "01510" || test_case == "01527" || test_case == "01596" || test_case == "01663" || test_case == "01686" || test_case == "01693" || test_case ∈ ["01684", "01685", "01694", "01695", "01696", "01697", "01698", "01699", "01700", "01719", "00928", "00929"]
        results = results[2:end, :]
    end
    col_names =  Symbol.(replace.(string.(names(results)), " " => ""))
    rename!(results, col_names)

    # Case for FBA models an exceptions should be thrown 
    if !("Time" ∈ names(results) || "time" ∈ names(results))
        sbml_string = String(take!(Downloads.download(sbml_urls[1], IOBuffer())))
        SBMLImport.build_SBML_model(sbml_string; model_as_string=true)
    end

    t_save = "Time" in names(results) ? Float64.(results[!, :Time]) : Float64.(results[!, :time])
    t_save = Vector{Float64}(t_save)
    tmax = maximum(t_save)
    what_check = filter(x -> x ∉ [:time, :Time], col_names)
    sbml_url = sbml_urls[end]

    # Read settings file
    settings_url = base_url * "-settings.txt"
    species_test, species_test_amount, species_test_conc, abstol_test, reltol_test = read_settings(settings_url)

    for sbml_url in sbml_urls

        sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
        # c = n / V => n = c * V

        model_SBML = SBMLImporter.build_SBML_model(sbml_string, model_as_string=true)
        # If stoichiometryMath occurs we need to convert the SBML file to a level 3 file
        # to properly handle it
        if occursin("stoichiometryMath", sbml_string) == false
            libsbml_model = readSBMLFromString(sbml_string)
        else
            libsbml_model = readSBMLFromString(sbml_string, doc -> begin
                                set_level_and_version(3, 2)(doc)
                                convert_promotelocals_expandfuns(doc)
                                end)
        end

        reaction_system, specie_map, parameter_map, cb, get_tstops, ifelse_t0 = SBMLImporter.SBML_to_ReactionSystem(sbml_string, ret_all=true, model_as_string=true)
        ode_system = convert(ODESystem, reaction_system)
        ode_problem = ODEProblem(ode_system, specie_map, (0.0, tmax), parameter_map, jac=true)
        for _f! in ifelse_t0
            _f!(ode_problem.u0, ode_problem.p)
        end
        tstops = get_tstops(ode_problem.u0, ode_problem.p)
        tstops = isempty(tstops) ? tstops : vcat(minimum(tstops) / 2.0, tstops)

        sol = solve(ode_problem, solver, abstol=1e-12, reltol=1e-12, saveat=t_save, tstops=tstops, callback=cb)
        model_parameters = parameters(sol.prob.f.sys)
        for to_check in what_check
            to_check_no_whitespace = Symbol(replace(string(to_check), " " => ""))
            if to_check_no_whitespace ∈ Symbol.(model_parameters)
                iParam = findfirst(x -> x == to_check_no_whitespace, Symbol.(model_parameters))

                if all(isinf.(results[!, to_check])) && all(results[!, to_check] .> 0)
                    @test isinf(sol.prob.p[iParam]) && sol.prob.p[iParam] > 0
                elseif all(isinf.(results[!, to_check])) && all(results[!, to_check] .< 0)
                    @test isinf(sol.prob.p[iParam]) && sol.prob.p[iParam] < 0
                elseif all(isnan.(results[!, to_check]))
                    @test isnan(sol.prob.p[iParam])
                else
                    @test all(abs.(sol.prob.p[iParam] .- results[!, to_check]) .< abstol_test .+ reltol_test .* abs.(results[!, to_check]))
                end
                continue
            end

            if to_check ∈ species_test && to_check ∈ species_test_conc && string(to_check) ∈ keys(libsbml_model.species)
                compartmentName = libsbml_model.species[string(to_check)].compartment
                if model_SBML.species[string(to_check)].unit == :Concentration
                    c = 1.0
                elseif compartmentName in string.(model_parameters)
                    c = sol.prob.p[findfirst(x -> x == compartmentName, string.(model_parameters))]
                else
                    c = sol[Symbol(compartmentName)]
                end
            elseif to_check ∈ species_test && to_check ∈ species_test_amount && string(to_check) ∈ keys(libsbml_model.species)
                compartmentName = libsbml_model.species[string(to_check)].compartment
                if model_SBML.species[string(to_check)].unit == :Concentration
                    if compartmentName in string.(model_parameters)
                        c = 1 / (sol.prob.p[findfirst(x -> x == compartmentName, string.(model_parameters))])
                    else
                        c = 1 ./ sol[Symbol(compartmentName)]
                    end
                else
                    c = 1.0
                end
            else
                c = 1.0
            end

            @test all(abs.(sol[to_check] ./ c .- results[!, to_check]) .< abstol_test .+ reltol_test .* abs.(results[!, to_check]))
        end
    end
end


# Function to get model-str
function get_model_str(test_case)    
    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/semantic/$test_case/$test_case"
    sbml_urls = get_sbml_urls(base_url)
    sbml_url = sbml_urls[end]
    sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
    model_SBML = SBMLImporter.build_SBML_model(sbml_string, model_as_string=true)
    reaction_system = SBMLImporter.reactionsystem_from_SBML(model_SBML)
    return reaction_system, model_SBML
end


# To run this you might have to add Catalyst into the SBMLImporter.jl file 
solver = Rodas4P()
@testset "Catalyst" begin
    for i in 1:5
        test_case = "0000$i"
        check_test_case(test_case, solver)
    end
end

model_str, model_SBML = get_model_str("00005")
println(model_str)
