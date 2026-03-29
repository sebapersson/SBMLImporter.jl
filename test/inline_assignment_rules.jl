using Catalyst, Downloads, ModelingToolkitBase, OrdinaryDiffEqRosenbrock,
    SBMLImporter, Test

# Boehm published model commonly used in benchmarks
# Inline assignment rules
path_SBML = joinpath(@__DIR__, "Models", "model_Boehm_JProteomeRes2014.xml")
rn1, cb1 = load_SBML(path_SBML; inline_assignment_rules = true)
u01, ps1 = get_u0_map(rn1), get_parameter_map(rn1)
oprob1 = ODEProblem(rn1, u01, (0.0, 10.0), ps1)
sol1 = solve(oprob1, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, saveat = 1:10)
# Use mtkcompile to inline assignment rules
rn2, cb2 = load_SBML(path_SBML; inline_assignment_rules = false)
u02, ps2 = get_u0_map(rn2), get_parameter_map(rn2)
sys = mtkcompile(ode_model(rn2))
oprob2 = ODEProblem(sys, merge(Dict(u02), Dict(ps2)), (0.0, 10.0))
sol2 = solve(oprob2, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, saveat = 1:10)
# Test equivalent solutions regardless of inline mode
for name in unknowns(sys)
    @test all(.≈(sol1[name], sol2[name], atol = 1.0e-6))
end

# Published model with lots of edge cases. The model has a log in it making it not possible
# to solve it, but we can test number of final unknowns
path_sbml = joinpath(@__DIR__, "Models", "BIOMD0000000627_url.xml")
rn, cb = load_SBML(
    path_sbml; inline_assignment_rules = true, ifelse_to_callback = true,
)
u0, ps = get_u0_map(rn), get_parameter_map(rn)
oprob = ODEProblem(rn, u0, (0.0, 10.0), ps)
@test length(oprob.u0) == 66

# Another published model with lots of edge cases.
# Inline assignment rules
path_SBML = joinpath(@__DIR__, "Models", "EDES_1_0.xml")
rn1, cb1 = load_SBML(path_SBML; inline_assignment_rules = true, ifelse_to_callback = false)
u01, ps1 = get_u0_map(rn1), get_parameter_map(rn1)
oprob1 = ODEProblem(rn1, u01, (0.0, 10.0), ps1)
sol1 = solve(oprob1, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, saveat = 1:10)
# Inline via mtkcompile
rn2, cb2 = load_SBML(path_SBML; inline_assignment_rules = false, ifelse_to_callback = false)
u02, ps2 = get_u0_map(rn2), get_parameter_map(rn2)
sys = mtkcompile(ode_model(rn2))
#u02 = first.(u02) .=> 0.0
oprob2 = ODEProblem(sys, merge(Dict(u02), Dict(ps2)), (0.0, 10.0))
sol2 = solve(oprob2, Rodas5P(), abstol = 1.0e-3, reltol = 1.0e-8, saveat = 1:10)
for name in unknowns(sys)
    @test all(.≈(sol1[name], sol2[name], atol = 1.0e-6))
end
