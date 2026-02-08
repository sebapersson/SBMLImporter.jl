using Catalyst, ModelingToolkitBase, OrdinaryDiffEqRosenbrock, SBMLImporter, Test

# With callback
path_SBML = joinpath(@__DIR__, "Models", "model_Brannmark_JBC2010.xml")
prn1, cb1 = load_SBML(path_SBML)
sys1 = mtkcompile(ode_model(prn1.rn))
oprob1 = ODEProblem(sys1, merge(Dict(prn1.u0), Dict(prn1.p)), (0.0, 10.0))
sol1 = solve(oprob1, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, callback = cb1)
# Without callback
prn2, cb2 = load_SBML(path_SBML; ifelse_to_callback = false)
sys2 = mtkcompile(ode_model(prn2.rn))
oprob2 = ODEProblem(sys2, merge(Dict(prn2.u0), Dict(prn2.p)), (0.0, 10.0))
sol2 = solve(oprob2, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8)
for name in unknowns(sys1)
    @test all(.≈(sol1[name], sol2[name], atol = 1.0e-9))
end
