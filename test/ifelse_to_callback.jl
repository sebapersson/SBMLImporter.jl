using SBMLImporter, ModelingToolkit, OrdinaryDiffEq, Test

# With callback
path_SBML = joinpath(@__DIR__, "Models", "model_Brannmark_JBC2010.xml")
prn1, cb1 = load_SBML(path_SBML)
sys1 = structural_simplify(convert(ODESystem, prn1.rn))
oprob1 = ODEProblem(sys1, prn1.u0, (0.0, 10.0), prn1.p)
sol1 = solve(oprob1, Rodas5P(), abstol=1e-3, reltol=1e-8, callback = cb1)
# Without callback
prn2, cb2 = load_SBML(path_SBML; ifelse_to_callback = false)
sys2 = structural_simplify(convert(ODESystem, prn2.rn))
oprob2 = ODEProblem(sys2, prn2.u0, (0.0, 10.0), prn2.p)
sol2 = solve(oprob2, Rodas5P(), abstol=1e-3, reltol=1e-8)
for name in states(sys1)
    all(.≈(sol1[name], sol2[name], atol = 1e-9))
end
