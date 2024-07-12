using OrdinaryDiffEq, SBMLImporter, ModelingToolkit, Catalyst, Test

path_model = joinpath(@__DIR__, "Models", "SBML", "Model.jl")
isfile(path_model) && rm(path_model)

path_SBML = joinpath(@__DIR__, "Models", "model_Boehm_JProteomeRes2014.xml")
parsed_rn, cb = load_SBML(path_SBML; write_to_file = true)
sys_ref = structural_simplify(convert(ODESystem, parsed_rn.rn))
oprob_ref = ODEProblem(sys_ref, parsed_rn.u₀, (0.0, 5.0), parsed_rn.p, jac = true)
sol_ref = solve(oprob_ref, Rodas4P(), saveat = 1:5)

include(path_model)

rn, u0map, pmap = get_reaction_system([])
sys_check = structural_simplify(convert(ODESystem, rn))
oprob_check = ODEProblem(sys_ref, u0map, (0.0, 5.0), pmap, jac = true)
sol_check = solve(oprob_check, Rodas4P(), saveat = 1:5)

@test oprob_check.u0 == oprob_ref.u0
@test oprob_check.p == oprob_ref.p
@test sol_ref.u == sol_check.u

rm(path_model)
rm(joinpath(@__DIR__, "Models", "SBML"))
