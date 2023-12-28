using SafeTestsets

@safetestset "SBML test-suite" begin
    include("Testsuite_catalyst.jl")
end