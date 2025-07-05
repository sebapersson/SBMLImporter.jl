#=
    Compare Catalyst vs ODEProblem import. The latter is needed by PEtab SciML, and
    rewrites a ModelSBML to an ODEProblem instead of rewriting to a ReactionSystem
=#
using SBMLImporter, ComponentArrays, OrdinaryDiffEqRosenbrock, Test

function test_odeproblem_import(path_SBML; atol = 1e-6)::Nothing
    model_SBML = SBMLImporter.parse_SBML(path_SBML, false; model_as_string = false,
                                         inline_assignment_rules=true)
    model_SBML_prob = SBMLImporter.ModelSBMLProb(model_SBML)
    @info "Model $(model_SBML.name)"

    prn1, cb1 = load_SBML(path_SBML; inline_assignment_rules = true)
    oprob1 = ODEProblem(prn1.rn, prn1.u0, (0.0, 10.0), prn1.p)
    sol1 = solve(oprob1, Rodas5P(), abstol = 1e-10, reltol = 1e-10, saveat = 1:10,
                 callback = cb1)
    oprob2, compute_u0, cb2 = SBMLImporter._get_odeproblem(model_SBML_prob, model_SBML)
    oprob2.u0 .= compute_u0(oprob2.p, 0.0)
    sol2 = solve(oprob2, Rodas5P(), abstol = 1e-10, reltol = 1e-10, saveat = 1:10,
                 callback = cb2)
    for (i, uid) in pairs(Symbol.(model_SBML_prob.umodel))
        @test all(.≈(sol2[i, :], sol1[uid], atol = atol))
    end
    return nothing
end

models = ["brusselator_kinetic_parameters.xml", "model_Boehm_JProteomeRes2014.xml",
          "model_Brannmark_JBC2010.xml", "model_Bachmann_MSB2011.xml",
          "model_Beer_MolBioSystems2014.xml", "model_Isensee_JCB2018.xml",
          "model_Fujita_SciSignal2010.xml", "model_Smith_BMCSystBiol2013.xml"]
@testset "ODEProblem import" begin
    for model in models
        test_odeproblem_import(joinpath(@__DIR__, "Models", model))
    end
end
