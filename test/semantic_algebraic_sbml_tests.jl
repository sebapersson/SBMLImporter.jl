using Catalyst, CSV, DataFrames, Downloads, ModelingToolkitBase, OrdinaryDiffEqBDF,
    OrdinaryDiffEqRosenbrock, SBML, SBMLImporter, Test

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "testsuite_support.jl"))

#=
    NOTE ON LICENSING

    ModelingToolkit v11 is expected to be (or may already be) AGPL-licensed, while
    ModelingToolkitBase remains MIT-licensed (see link below).

    SBMLImporter uses ModelingToolkit **only** in the test suite, to verify that an
    imported `ODESystem` can be converted into a form that is simulatable when a user
    chooses to work with AGPL-licensed ModelingToolkit functionality. This only applies
    to a small subset of SBML tests with algebraic rules. ModelingToolkit is also tested
    separately to ensure nothing breaks in the simplification using only
    ModelingToolkitBase

    SBMLImporter does not ship ModelingToolkit functionality at runtime; the package
    itself depends only on ModelingToolkitBase.

    See discussion on ModelingToolkit getting APGL-license:
    https://discourse.julialang.org/t/modelingtoolkit-v11-library-split-and-licensing-community-feedback-requested/134396
=#
import ModelingToolkit

solver = Rodas4P()
@testset "Semantic algebraic test-suite" begin
    for test_case in ALGEBRAIC_TESTS
        fbdf_tests = [
            "00028", "00173", "00269", "01260", "01504", "01505", "01506", "01669",
            "01670", "01671", "00946",
        ]
        if test_case in fbdf_tests
            check_semantic_test_case(test_case, FBDF())
            continue
        end
        check_semantic_test_case(test_case, solver)
    end
end
