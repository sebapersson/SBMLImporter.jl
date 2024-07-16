using SBMLImporter, ModelingToolkit, OrdinaryDiffEq, Downloads, Test

# Boehm published benchmark model
path_SBML = joinpath(@__DIR__, "Models", "model_Boehm_JProteomeRes2014.xml")
prn1, cb1 = load_SBML(path_SBML; inline_assignment_rules = true)
oprob1 = ODEProblem(prn1.rn, prn1.u0, (0.0, 10.0), prn1.p)
sol1 = solve(oprob1, Rodas5P(), abstol=1e-3, reltol=1e-8)
prn2, cb2 = load_SBML(path_SBML; inline_assignment_rules = false)
sys = structural_simplify(convert(ODESystem, prn2.rn))
oprob2 = ODEProblem(sys, prn2.u0, (0.0, 10.0), prn2.p)
sol2 = solve(oprob2, Rodas5P(), abstol=1e-3, reltol=1e-8)
for name in states(sys)
    all(.≈(sol1[name], sol2[name], atol = 1e-9))
end

# Published model with lots of edge cases. The model has a log in it making it not possible
# to solve it, but we can test number of final states
sbml_url = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000627.3?filename=BIOMD0000000627_url.xml"
sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
mdl, cb = load_SBML(sbml_string; inline_assignment_rules = true, ifelse_to_callback = true,
                    model_as_string = true)
oprob = ODEProblem(mdl.rn, mdl.u0, (0.0, 10.0), mdl.p)
@test length(oprob.u0) == 66

# Another published model with lots of edge cases.
path_SBML = joinpath(@__DIR__, "Models", "EDES_1_0.xml")
prn1, cb1 = load_SBML(path_SBML; inline_assignment_rules = true, ifelse_to_callback = false)
oprob1 = ODEProblem(prn1.rn, prn1.u0, (0.0, 10.0), prn1.p)
sol1 = solve(oprob1, Rodas5P(), abstol=1e-3, reltol=1e-8, saveat=1:10)
prn2, cb2 = load_SBML(path_SBML; inline_assignment_rules = false, ifelse_to_callback = false)
sys = structural_simplify(convert(ODESystem, prn2.rn))
oprob2 = ODEProblem(sys, prn2.u0, (0.0, 10.0), prn2.p)
sol2 = solve(oprob2, Rodas5P(), abstol=1e-3, reltol=1e-8, saveat=1:10)
for name in states(sys)
    all(.≈(sol1[name], sol2[name], atol = 1e-9))
end
