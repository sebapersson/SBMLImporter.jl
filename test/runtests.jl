using SafeTestsets

core_only = get(ENV, "CORE_ONLY", "false") == "true"

@safetestset "Aqua Quality Check" begin
    include("aqua.jl")
end

@safetestset "Currated jump simulation test" begin
    include("curated_jump_test.jl")
end

@safetestset "Inline assignment rules" begin
    include("inline_assignment_rules.jl")
end

@safetestset "Kinetic law parameters" begin
    include("kineticlaw_parameters.jl")
end

@safetestset "Ifelse to callback" begin
    include("ifelse_to_callback.jl")
end

@safetestset "Test write parsed model to file" begin
    include("write_model.jl")
end

@safetestset "Reaction metadata" begin
    include("metadata.jl")
end

if !core_only
    @safetestset "Import as ODEProblem" begin
        include("odeproblem_import.jl")
    end

    @safetestset "SBML large models" begin
        include("large_models.jl")
    end

    @safetestset "SBML semantic test-suite" begin
        include("semantic_sbml_tests.jl")
    end

    @safetestset "SBML stochastic test-suite" begin
        include("stochastic_sbml_tests.jl")
    end
end
