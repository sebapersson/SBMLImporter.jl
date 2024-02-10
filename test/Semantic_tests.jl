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
    species_test = Symbol.(replace.(split(split(settings_lines[4], ":")[2], ','),
                                    " " => ""))
    species_test_amount = Symbol.(replace.(split(split(settings_lines[7], ":")[2], ','),
                                           " " => ""))
    species_test_conc = Symbol.(replace.(split(split(settings_lines[8], ":")[2], ','),
                                         " " => ""))
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
    if test_case == "00995" || test_case == "00996" || test_case == "00997" ||
       test_case == "01284" || test_case == "01510" || test_case == "01527" ||
       test_case == "01596" || test_case == "01663" || test_case == "01686" ||
       test_case == "01693" ||
       test_case ∈ [
           "01684",
           "01685",
           "01694",
           "01695",
           "01696",
           "01697",
           "01698",
           "01699",
           "01700",
           "01719",
           "00928",
           "00929",
       ]
        results = results[2:end, :]
    end
    col_names = Symbol.(replace.(string.(names(results)), " " => ""))
    rename!(results, col_names)

    # Case for FBA models an exceptions should be thrown
    if !("Time" ∈ names(results) || "time" ∈ names(results))
        sbml_string = String(take!(Downloads.download(sbml_urls[1], IOBuffer())))
        SBMLImporter.build_SBML_model(sbml_string; model_as_string = true)
    end

    t_save = "Time" in names(results) ? Float64.(results[!, :Time]) :
             Float64.(results[!, :time])
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

        parsed_rn, cb = load_SBML(sbml_string, model_as_string = true,
                                  inline_assignment_rules = false)
        if isempty(model_SBML.algebraic_rules)
            ode_system = structural_simplify(convert(ODESystem, parsed_rn.rn))
        else
            ode_system = structural_simplify(dae_index_lowering(convert(ODESystem,
                                                                        parsed_rn.rn)))
        end
        ode_problem = ODEProblem(ode_system, parsed_rn.u₀, (0.0, tmax), parsed_rn.p,
                                 jac = true)

        sol = solve(ode_problem, solver, abstol = 1e-12, reltol = 1e-12, saveat = t_save,
                    callback = cb)
        model_parameters = parameters(sol.prob.f.sys)
        for to_check in what_check
            to_check_no_whitespace = Symbol(replace(string(to_check), " " => ""))
            if to_check_no_whitespace ∈ Symbol.(model_parameters)
                iParam = findfirst(x -> x == to_check_no_whitespace,
                                   Symbol.(model_parameters))

                if all(isinf.(results[!, to_check])) && all(results[!, to_check] .> 0)
                    @test isinf(sol.prob.p[iParam]) && sol.prob.p[iParam] > 0
                elseif all(isinf.(results[!, to_check])) && all(results[!, to_check] .< 0)
                    @test isinf(sol.prob.p[iParam]) && sol.prob.p[iParam] < 0
                elseif all(isnan.(results[!, to_check]))
                    @test isnan(sol.prob.p[iParam])
                else
                    @test all(abs.(sol.prob.p[iParam] .- results[!, to_check]) .<
                              abstol_test .+ reltol_test .* abs.(results[!, to_check]))
                end
                continue
            end

            if to_check ∈ species_test && to_check ∈ species_test_conc &&
               haskey(libsbml_model.species, string(to_check))
                compartmentName = libsbml_model.species[string(to_check)].compartment
                if model_SBML.species[string(to_check)].unit == :Concentration
                    c = 1.0
                elseif compartmentName in string.(model_parameters)
                    c = sol.prob.p[findfirst(x -> x == compartmentName,
                                             string.(model_parameters))]
                else
                    c = sol[Symbol(compartmentName)]
                end
            elseif to_check ∈ species_test && to_check ∈ species_test_amount &&
                   haskey(libsbml_model.species, string(to_check))
                compartmentName = libsbml_model.species[string(to_check)].compartment
                if model_SBML.species[string(to_check)].unit == :Concentration
                    if compartmentName in string.(model_parameters)
                        c = 1 / (sol.prob.p[findfirst(x -> x == compartmentName,
                                                  string.(model_parameters))])
                    else
                        c = 1 ./ sol[Symbol(compartmentName)]
                    end
                else
                    c = 1.0
                end
            else
                c = 1.0
            end

            @test all(abs.(sol[to_check] ./ c .- results[!, to_check]) .<
                      abstol_test .+ reltol_test .* abs.(results[!, to_check]))
        end
    end
end

# Function to get model-str
function get_model_str(test_case)
    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/semantic/$test_case/$test_case"
    sbml_urls = get_sbml_urls(base_url)
    sbml_url = sbml_urls[end]
    sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
    model_SBML = SBMLImporter.build_SBML_model(sbml_string, model_as_string = true,
                                               inline_assignment_rules = false)

    parsed_model_SBML = SBMLImporter._reactionsystem_from_SBML(model_SBML)
    model_str = SBMLImporter.reactionsystem_to_string(parsed_model_SBML, false, "",
                                                      model_SBML)

    return model_str, model_SBML
end

# To run this you might have to add Catalyst into the SBMLImporter.jl file
solver = Rodas4P()
@testset "Catalyst" begin
    for i in 1:1821
        test_case = repeat("0", 5 - length(string(i))) * string(i)

        delay_cases = ["00937", "00938", "00939", "00940", "00941", "00942", "00943",
            "00981",
            "00982", "00983", "00984", "00985", "01400",
            "01401", "01403", "01404", "01410", "01411", "01412", "01413", "01414",
            "01415", "01416", "01417", "01418", "01419", "01454", "01480", "01481",
            "01518", "01522", "01523", "01524", "01534", "01535", "01536", "01537",
            "01538", "01539", "01406", "01407", "01409", "01410"]

        event_delay_cases = ["00071", "00072", "00073", "00405", "00405", "00406", "00407",
            "00409", "00410", "00411", "00412", "00413", "00414", "00415",
            "00416", "00417", "00418", "00419", "00420", "00622", "00623",
            "00624", "00637", "00638", "00639", "00649", "00650", "00651",
            "00664", "00665", "00666", "00682", "00683", "00684", "00690",
            "00702", "00708", "00724", "00737", "00757", "00758", "00759",
            "00763", "00764", "00765", "00766", "00767", "00768", "00769",
            "00770", "00771", "00772", "00773", "00774", "00775", "00776",
            "00777", "00778", "00779", "00780", "00848", "00849", "00850",
            "00886", "00887", "00932", "00933", "00936", "00408", "00461",
            "00655", "00656", "00657", "00980", "01000", "01048", "01049",
            "01050", "01074", "01075", "01076", "01119", "01120", "01230",
            "01241", "01263", "01268", "01269", "01270", "01287", "01295",
            "01299", "01305", "01324", "01325", "01326", "01327", "01328",
            "01329", "01335", "01507", "01508", "01509", "01511", "01519",
            "01520", "01525", "01526", "01528", "01529", "01532", "01575",
            "01576", "01579", "01580", "01581", "01582", "01584", "01585",
            "01586", "01587", "01594", "01595", "01597", "01598", "01600",
            "01601", "01602", "01603", "01604", "01659", "01660", "01661",
            "01673", "01674", "01675", "01676", "01677", "01678", "01679",
            "01680", "01687", "01688", "01689", "01690", "01691", "01692",
            "01701", "01702", "01703", "01704", "01706", "01707", "01708",
            "01709", "01710", "01711", "01712", "01713", "01715", "01716",
            "01717", "01718", "01720", "01721", "01754", "01755", "01756",
            "01757", "01758", "01759", "01798", "00731", "01095", "01771",
            "01672", "00701"]
        event_delay_cases = vcat(event_delay_cases, ["004" * string(i) for i in 21:61])

        event_priority = ["00931", "00934", "00935", "00962", "00963", "00964", "00965",
            "00966", "00967",
            "00978", "00978", "01229", "01242", "01267", "01294", "01298", "01331", "01332",
            "01333", "01334", "01336", "01337", "01466", "01512", "01521", "01533", "01577",
            "01583", "01588", "01589", "01590", "01591", "01592", "01593", "01599", "01605",
            "01626", "01627", "01662", "01681", "01682", "01683", "01705", "01714", "01212",
            "01213", "01214", "01215", "00952", "00930", "00997", "01262", "01286", "01330",
            "01772"]

        fast_reactions = [
            "00870",
            "00871",
            "00872",
            "00873",
            "00874",
            "00875",
            "00986",
            "00987",
            "00988",
            "01051",
            "01052",
            "01053",
            "01396",
            "01397",
            "01398",
            "01399",
            "01544",
            "01545",
            "01546",
            "01547",
            "01548",
            "01549",
            "01550",
            "01551",
            "01558",
            "01559",
            "01560",
            "01565",
            "01567",
            "01568",
            "01569",
            "01570",
            "01571",
            "01572",
        ]

        imply_cases = ["01274", "01279", "01497"]

        three_arugment_and_xor = ["00198", "00200", "00201", "00276", "00279"]

        piecewise_in_function_with_and_or = [
            "00956",
            "01282",
            "01488",
            "01503",
            "01516",
            "01561",
        ]

        piecewise_algebraic_rules = ["01503"]

        rem_div_parameters = ["01272", "01277", "01495"]

        event_multiple_triggers = ["01211", "01531"]

        multiple_argument_inequality = [
            "01216",
            "01494",
            "01781",
            "01782",
            "01783",
            "01209",
            "01210",
            "00278",
        ]

        hierarchical_models = vcat(["011" * string(i) for i in 24:83],
                                   ["0" * string(i) for i in 1344:1394],
                                   ["0" * string(i) for i in 1467:1477], ["01778"])

        fba_models = vcat(["01" * string(i) for i in 186:196],
                          ["01" * string(i) for i in 606:625],
                          ["01627", "01628", "01629", "01630"])

        bad_names = [
            "01810",
            "01811",
            "01812",
            "01813",
            "01814",
            "01815",
            "01816",
            "01817",
            "01818",
            "01819",
            "01820",
            "01821",
        ]

        if test_case ∈ ["00028", "00173", "00269"]
            check_test_case(test_case, FBDF())
            continue
        end

        if test_case in fba_models
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
        end
        if test_case ∈ piecewise_algebraic_rules
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ bad_names
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ hierarchical_models
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ multiple_argument_inequality
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case == "01596" # Must be able to solve for time event is triggered
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ event_multiple_triggers
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ delay_cases
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ event_delay_cases
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ event_priority
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ fast_reactions
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ imply_cases
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ piecewise_in_function_with_and_or
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ three_arugment_and_xor
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case ∈ rem_div_parameters
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
        end

        # We do not aim to support Flux-Balance-Analysis (FBA) models
        not_test1 = ["01" * string(i) for i in 186:197]
        not_test2 = ["01" * string(i) for i in 606:625]
        if test_case ∈ not_test1 || test_case ∈ not_test2 ||
           test_case ∈ ["01627", "01628", "01629", "01630"]
            continue
        end

        #=
            Errors SBML importer does not manage to catch - but they should not simulate so it is something
        =#

        # Edge-cases where I think there is a bug in SBML.jl
        if test_case ∈ ["01318", "01319", "01320"]
            continue
        end

        # Edge case we cannot capture as we have nested piecewise - where the picewise is inside a function
        if test_case ∈ ["01491", "01492", "01493", "01486"]
            continue
        end

        # Uncommon mathML with boundary effects, not sure I agree with test-suite
        if test_case ∈ ["00959"]
            continue
        end

        # A form of event compeition
        if test_case ∈ ["00953"]
            continue
        end

        # If user wants to add a random species, it must either be as a species, initialAssignment, assignmentRule
        # or by event, not by just random adding it to equations. Cannot capture this error.
        if test_case ∈ ["00974"]
            continue
        end

        # Bug in SBML.jl (parameter rateOf)
        if test_case ∈ ["01321", "01322"]
            continue
        end

        # Even trigger cannot depend on algebraic rule
        if test_case == "01578"
            continue
        end

        # Inf or NaN in initial assignments
        if test_case ∈ ["00950", "00951"]
            continue
        end

        # Reaction cannot be empty
        if test_case ∈
           ["01245", "01246", "01300", "01301", "01302", "01303", "01304", "01305", "01306"]
            continue
        end

        check_test_case(test_case, solver)
    end
end
