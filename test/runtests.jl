using SafeTestsets

@safetestset "SBML large model" begin
    include("Large_model.jl")
end

@safetestset "SBML semantic test-suite" begin
    include("Testsuite_catalyst.jl")
end

@safetestset "SBML stochastic test-suite" begin
    include("Stochastic_tests.jl")
end
