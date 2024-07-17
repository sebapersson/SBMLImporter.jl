using SafeTestsets

@safetestset "Aqua Quality Check" begin
    include("aqua.jl")
end

@safetestset "Currated jump simulation test" begin
    include("curated_jump_test.jl")
end

@safetestset "Inline assignment rules" begin
    include("inline_assignment_rules.jl")
end

@safetestset "Ifelse to callback" begin
    include("ifelse_to_callback.jl")
end

@safetestset "SBML large models" begin
    include("large_models.jl")
end

@safetestset "Test write parsed model to file" begin
    include("write_model.jl")
end

@safetestset "SBML semantic test-suite p1" begin
    include("semantic_sbml_tests_p1.jl")
end

@safetestset "SBML semantic test-suite p2" begin
    include("semantic_sbml_tests_p2.jl")
end

@safetestset "SBML semantic test-suite p3" begin
    include("semantic_sbml_tests_p3.jl")
end

@safetestset "SBML semantic test-suite p4" begin
    include("semantic_sbml_tests_p4.jl")
end

@safetestset "SBML stochastic test-suite" begin
    include("stochastic_sbml_tests.jl")
end
