using Catalyst, ModelingToolkitBase, OrdinaryDiffEqRosenbrock, SBMLImporter, Test

path_model = joinpath(@__DIR__, "Models", "SBML", "Boehm_JProteomeRes2014.jl")
isfile(path_model) && rm(path_model)

path = joinpath(@__DIR__, "Models", "model_Boehm_JProteomeRes2014.xml")
rn, cb = load_SBML(path; write_to_file = true)
u0, ps = get_u0_map(rn), get_parameter_map(rn)
sys_ref = mtkcompile(ode_model(rn))
oprob_ref = ODEProblem(sys_ref, merge(Dict(u0), Dict(ps)), (0.0, 5.0), jac = true)
sol_ref = solve(oprob_ref, Rodas4P(), saveat = 1:5)

include(path_model)

rn, u0map, pmap = get_reaction_system([])
rn = Catalyst.complete(rn)
sys_check = mtkcompile(ode_model(rn))
oprob_check = ODEProblem(sys_ref, merge(Dict(u0map), Dict(pmap)), (0.0, 5.0), jac = true)
sol_check = solve(oprob_check, Rodas4P(), saveat = 1:5)

@test oprob_check.u0 == oprob_ref.u0
@test oprob_check.p == oprob_ref.p
@test sol_ref.u == sol_check.u

rm(path_model)
rm(joinpath(@__DIR__, "Models", "SBML"))
