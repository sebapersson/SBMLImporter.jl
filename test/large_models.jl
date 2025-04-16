#=
    Test that import time for a larger model of 1265 species and 1000+ reactions
=#

using SBMLImporter, JumpProcesses, Downloads, ModelingToolkit, Test

path = joinpath(@__DIR__, "Models", "large_model1.xml")
b1 = @elapsed prn, cb = load_SBML(path; inline_assignment_rules = true)
# Usually takes around 4s locally, but better to brace for GitHub CI
@test b1 ≤ 40

# Large fceri-γ model with ≈ 58,000 reactions. Should take around 4min locally, but it
# is always good to brace for GitHub CI
path = joinpath(@__DIR__, "Models", "fceri_gamma2.xml")
b2 = @elapsed prn, cb = load_SBML(path; inline_assignment_rules = true)
@test b2 ≤ 500

# Test that SBMLImporter correctly enforces mass action
path = joinpath(@__DIR__, "Models", "egfr_net.xml")
b3 = @elapsed model, cb = load_SBML(path; massaction = true, inline_assignment_rules = true)
dprob = DiscreteProblem(model.rn, model.u0, (0.0, 0.0), model.p)
dprob = remake(dprob, u0 = Int64.(dprob.u0));
jprob = JumpProblem(model.rn, dprob, RSSA(), save_positions = (false, false))
@test b3 ≤ 20
@test length(jprob.massaction_jump.net_stoch) == 3749
