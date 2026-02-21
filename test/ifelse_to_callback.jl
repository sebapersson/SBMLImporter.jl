using Catalyst, ModelingToolkitBase, OrdinaryDiffEqRosenbrock, SBMLImporter, Test

# With callback
path_SBML = joinpath(@__DIR__, "Models", "model_Brannmark_JBC2010.xml")
rn1, cb1 = load_SBML(path_SBML)
u01, ps1 = get_u0_map(rn1), get_parameter_map(rn1)
sys1 = mtkcompile(ode_model(rn1))
oprob1 = ODEProblem(sys1, merge(Dict(u01), Dict(ps1)), (0.0, 10.0))
sol1 = solve(oprob1, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, callback = cb1)

# Without callback
rn2, cb2 = load_SBML(path_SBML; ifelse_to_callback = false)
u02, ps2 = get_u0_map(rn2), get_parameter_map(rn2)
sys2 = mtkcompile(ode_model(rn2))
oprob2 = ODEProblem(sys2, merge(Dict(u02), Dict(ps2)), (0.0, 10.0))
sol2 = solve(oprob2, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8)
for name in unknowns(sys1)
    @test all(.≈(sol1[name], sol2[name], atol = 1.0e-9))
end
