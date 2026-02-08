using Catalyst, Downloads, ModelingToolkitBase, OrdinaryDiffEqRosenbrock,
    SBMLImporter, Test

# Boehm published model commonly used in benchmarks
# Inline assignment rules
path_SBML = joinpath(@__DIR__, "Models", "model_Boehm_JProteomeRes2014.xml")
prn1, cb1 = load_SBML(path_SBML; inline_assignment_rules = true)
oprob1 = ODEProblem(prn1.rn, prn1.u0, (0.0, 10.0), prn1.p)
sol1 = solve(oprob1, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, saveat = 1:10)
# Use mtkcompile to inline assignment rules
prn2, cb2 = load_SBML(path_SBML; inline_assignment_rules = false)
sys = mtkcompile(ode_model(prn2.rn))
oprob2 = ODEProblem(sys, merge(Dict(prn2.u0), Dict(prn2.p)), (0.0, 10.0))
sol2 = solve(oprob2, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, saveat = 1:10)
# Test equivalent solutions regardless of inline mode
for name in unknowns(sys)
    @test all(.≈(sol1[name], sol2[name], atol = 1.0e-6))
end

# Published model with lots of edge cases. The model has a log in it making it not possible
# to solve it, but we can test number of final unknowns
sbml_url = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000627.3?filename=BIOMD0000000627_url.xml"
sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
prn, cb = load_SBML(
    sbml_string; inline_assignment_rules = true, ifelse_to_callback = true,
    model_as_string = true
)
oprob = ODEProblem(prn.rn, prn.u0, (0.0, 10.0), prn.p)
@test length(oprob.u0) == 66

# Another published model with lots of edge cases.
# Inline assignment rules
path_SBML = joinpath(@__DIR__, "Models", "EDES_1_0.xml")
prn1, cb1 = load_SBML(path_SBML; inline_assignment_rules = true, ifelse_to_callback = false)
oprob1 = ODEProblem(prn1.rn, prn1.u0, (0.0, 10.0), prn1.p)
sol1 = solve(oprob1, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, saveat = 1:10)
# Inline via mtkcompile
prn2, cb2 = load_SBML(path_SBML; inline_assignment_rules = false, ifelse_to_callback = false)
sys = mtkcompile(ode_model(prn2.rn))
u0 = first.(prn2.u0) .=> 0.0
oprob2 = ODEProblem(sys, merge(Dict(prn2.u0), Dict(prn2.p)), (0.0, 10.0))
sol2 = solve(oprob2, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, saveat = 1:10)
for name in unknowns(sys)
    @test all(.≈(sol1[name], sol2[name], atol = 1.0e-6))
end
