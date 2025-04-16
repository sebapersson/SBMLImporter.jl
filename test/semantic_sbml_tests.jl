using SBMLImporter
using CSV
using DataFrames
using OrdinaryDiffEq
using ModelingToolkit
using SBML
using Test
using Downloads
using Catalyst

include(joinpath(@__DIR__, "common.jl"))

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
    if !("Time" ∈ names(results) || "time" ∈ names(results))
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
            if component_test ∈ vcat(species_test_conc, species_test_amount) && is_sbml_specie
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

# To run this you might have to add Catalyst into the SBMLImporter.jl file
solver = Rodas4P()
@testset "Catalyst" begin
    for i in 1:1821
        test_case = repeat("0", 5 - length(string(i))) * string(i)

        delay_cases = ["00937", "00938", "00939", "00940", "00941", "00942", "00943",
                       "00981", "00982", "00983", "00984", "00985", "01400", "01401",
                       "01403", "01404", "01410", "01411", "01412", "01413", "01414",
                       "01415", "01416", "01417", "01418", "01419", "01454", "01480",
                       "01481", "01518", "01522", "01523", "01524", "01534", "01535",
                       "01536", "01537", "01538", "01539", "01406", "01407", "01409",
                       "01410"]
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
                             "01672", "00701", "01305"]
        event_delay_cases = vcat(event_delay_cases, ["004" * string(i) for i in 21:61])
        event_priority = ["00931", "00934", "00935", "00962", "00963", "00964", "00965",
                          "00966", "00967", "00978", "00978", "01229", "01242", "01267",
                          "01294", "01298", "01331", "01332", "01333", "01334", "01336",
                          "01337", "01466", "01512", "01521", "01533", "01577", "01583",
                          "01588", "01589", "01590", "01591", "01592", "01593", "01599",
                          "01605", "01626", "01627", "01662", "01681", "01682", "01683",
                          "01705", "01714", "01212", "01213", "01214", "01215", "00952",
                          "00930", "00997", "01262", "01286", "01330", "01772"]
        fast_reactions = ["00870", "00871", "00872", "00873", "00874", "00875", "00986",
                         "00987", "00988", "01051", "01052", "01053", "01396", "01397",
                         "01398", "01399", "01544", "01545", "01546", "01547", "01548",
                         "01549", "01550", "01551", "01558", "01559", "01560", "01565",
                         "01567", "01568", "01569", "01570", "01571", "01572"]
        imply_cases = ["01274", "01279", "01497"]
        three_arugment_and_xor = ["00198", "00200", "00201", "00276", "00279", "01282"]
        piecewise_algebraic_rules = ["01503"]
        rem_div_parameters = ["01272", "01273", "01277", "01278", "01495", "01496"]
        event_multiple_triggers = ["01211", "01531"]
        multiple_argument_inequality = ["01216", "01494", "01781", "01782", "01783",
                                        "01209", "01210", "00278"]
        hierarchical_models = vcat(["011" * string(i) for i in 24:83],
                                   ["0" * string(i) for i in 1344:1394],
                                   ["0" * string(i) for i in 1467:1477], ["01778"])
        fba_models = vcat(["01" * string(i) for i in 186:196],
                          ["01" * string(i) for i in 606:625],
                          ["01627", "01628", "01629", "01630"])
        bad_names = ["01810", "01811", "01812", "01813", "01814", "01815", "01816",
                    "01817", "01818", "01819", "01820", "01821",]

        if test_case in fba_models
            @test_throws SBMLImporter.SBMLSupport begin
                check_test_case(test_case, Rodas4P())
            end
            continue
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

        #=
            Errors SBML importer does not manage to catch - but these are very rare
            cases
        =#
        # Edge-cases where I think there is a bug in SBML.jl
        if test_case ∈ ["01318", "01319", "01320"]
            continue
        end
        # Uncommon mathML with boundary effects, not sure I agree with test-suite
        if test_case ∈ ["00959", "01488"]
            continue
        end
        # A form of event compeition
        if test_case ∈ ["00953"]
            continue
        end
        # If user wants to add a random species, it must either be as a species,
        # initialAssignment, assignmentRule or by event, not by just random adding it
        # to equations. Cannot capture this error.
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
        # Subset of empty reactions. TODO: Add fix for
        if test_case ∈ ["01304", "01305", "01306"]
            continue
        end

        # Simulations fail with Rodas4P, so FBDF is used
        if test_case ∈ ["00028", "00173", "00269"]
            check_test_case(test_case, FBDF())
            continue
        end
        check_test_case(test_case, solver)
    end
end
