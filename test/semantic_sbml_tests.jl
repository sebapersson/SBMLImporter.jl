using Catalyst, CSV, DataFrames, Downloads, ModelingToolkitBase, OrdinaryDiffEqBDF,
    OrdinaryDiffEqRosenbrock, SBML, SBMLImporter, Test

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "testsuite_support.jl"))

# Get test cases currently not passing
cases_not_passing = get_semantic_not_pass()
keys_cases_catch = filter(x -> x != :not_captured, keys(cases_not_passing))
cases_not_passing_catch = reduce(vcat, [cases_not_passing[k] for k in keys_cases_catch])
solver = Rodas4P()
@testset "Semantic test-suite" begin
    for i in 1:1821
        test_case = repeat("0", 5 - length(string(i))) * string(i)

        test_case in ALGEBRAIC_TESTS && continue

        if test_case in cases_not_passing_catch
            @test_throws SBMLImporter.SBMLSupport begin
                check_semantic_test_case(test_case, Rodas4P())
            end
            continue
        end
        if test_case in cases_not_passing[:not_captured]
            continue
        end
        if test_case == "01596" # Must be able to solve for time event is triggered
            @test_throws SBMLImporter.SBMLSupport begin
                check_semantic_test_case(test_case, Rodas4P())
            end
            continue
        end

        # Simulations fail with Rodas4P, so FBDF is used
        fbdf_tests = [
            "00028", "00173", "00269", "01260", "01504", "01505", "01506", "01669",
            "01670", "01671", "00946"
        ]
        if test_case in fbdf_tests
            check_semantic_test_case(test_case, FBDF())
            continue
        end
        check_semantic_test_case(test_case, solver)
    end
end
