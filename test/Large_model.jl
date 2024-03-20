#=
    Test that import time for a larger model of 1265 species and 1000+ reactions
=#

using SBMLImporter
using JumpProcesses
using Downloads
using ModelingToolkit
using Test

sbml_url = "https://www.ebi.ac.uk/biomodels/model/download/MODEL1112100000.2?filename=MODEL1112100000_url.xml"
sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))

b1 = @elapsed parsed_rn, cb = load_SBML(sbml_string; model_as_string = true)

# Usually takes around 4s locally, but better to brace for GitHub CI
@test b1 ≤ 20

# Large fceri-γ model with ≈ 58,000 reactions. Should take around 4min locally, but it
# is always good to brace for GitHub CI
path_SBML = joinpath(@__DIR__, "Models", "fceri_gamma2.xml")
b2 = @elapsed parsed_rn, cb = load_SBML(path_SBML)
@test b2 ≤ 500

# Test that SBMLImporter correctly enforces mass action
path_SBML = joinpath(@__DIR__, "Models", "egfr_net.xml")
b3 = @elapsed model, cb = load_SBML(path_SBML; mass_action=true)
dprob = DiscreteProblem(model.rn, model.u₀, (0.0,0.0), model.p)
dprob = remake(dprob, u0 = Int64.(dprob.u0));
jprob = JumpProblem(model.rn, dprob, RSSA(), save_positions=(false,false))
@test b3 ≤ 20
@test length(jprob.massaction_jump.net_stoch) == 3749

# Lots of edge cases
sbml_url = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000627.3?filename=BIOMD0000000627_url.xml"
sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
mdl, cb = load_SBML(sbml_string; inline_assignment_rules=false, ifelse_to_callback=true, check_massaction=false,
                    model_as_string=true)
sys = structural_simplify(convert(ODESystem, mdl.rn))
@test length(states(sys)) == 66
