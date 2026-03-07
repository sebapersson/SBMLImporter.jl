using Catalyst, OrdinaryDiffEqRosenbrock, SBMLImporter, Test

#=
    Test that if kinetic law parameters are inlined correctly if the user selects that
    option. Further, tests different duplicate id cases.
=#

# First case with no duplicate ids
path_SBML = joinpath(@__DIR__, "Models", "brusselator_kinetic_parameters.xml")
rn1, cb1 = load_SBML(path_SBML; inline_kineticlaw_parameters = false)
u01, ps1 = get_u0_map(rn1), get_parameter_map(rn1)
@test string.(parameters(rn1)) == ["B", "A", "k1_r2", "k1_r4", "compartment"]
@test string(ps1[3]) == "k1_r2 => 1.0"
@test string(ps1[4]) == "k1_r4 => 1.0"

rn2, cb2 = load_SBML(path_SBML; inline_kineticlaw_parameters = true)
u02, ps2 = get_u0_map(rn2), get_parameter_map(rn2)
@test string.(parameters(rn2)) == ["B", "A", "compartment"]
oprob1 = ODEProblem(rn1, u01, (0.0, 10.0), ps1)
sol1 = solve(oprob1, Rodas5P())
oprob2 = ODEProblem(rn2, u02, (0.0, 10.0), ps2)
sol2 = solve(oprob2, Rodas5P())
@test sum(reduce(vcat, sol1.u .- sol2.u) .^ 2) ≤ 1.0e-10

# With duplicate id parameters, where the parameters have the same value and therefore
# the model can be imported
path_SBML = joinpath(@__DIR__, "Models", "brusselator.xml")
rn, cb = load_SBML(path_SBML; inline_kineticlaw_parameters = false)
@test string.(parameters(rn)) == ["B", "A", "k1", "C"]
@test string(get_parameter_map(rn)[3]) == "k1 => 1.0"

# With duplicate id parameters, where the parameters have the same value and therefore
# the model cannot be imported
path_SBML = joinpath(@__DIR__, "Models", "brusselator_kinetic_parameters_fail.xml")
@test_throws SBMLImporter.SBMLSupport begin
    _, _ = load_SBML(path_SBML; inline_kineticlaw_parameters = false)
end
