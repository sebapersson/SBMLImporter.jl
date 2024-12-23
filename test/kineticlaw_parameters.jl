using SBMLImporter, Catalyst, OrdinaryDiffEq, Test

#=
    Test that if kinetic law parameters are inlined correctly if the user selects that
    option. Further, tests different duplicate id cases.
=#

# First case with no duplicate ids
path_SBML = joinpath(@__DIR__, "Models", "brusselator_kinetic_parameters.xml")
prn1, cb1 = load_SBML(path_SBML; inline_kineticlaw_parameters = false)
@test string.(parameters(prn1.rn)) == ["B", "A", "k1_r2", "k1_r4", "C"]
@test string(prn1.p[3]) == "k1_r2 => 1.0"
@test string(prn1.p[4]) == "k1_r4 => 1.0"
prn2, cb2 = load_SBML(path_SBML; inline_kineticlaw_parameters = true)
@test string.(parameters(prn2.rn)) == ["B", "A", "C"]
oprob1 = ODEProblem(prn1.rn, prn1.u0, (0.0, 10.0), prn1.p)
sol1 = solve(oprob1, Rodas5P())
oprob2 = ODEProblem(prn2.rn, prn2.u0, (0.0, 10.0), prn2.p)
sol2 = solve(oprob2, Rodas5P())
@test sum(reduce(vcat, sol1.u .- sol2.u).^2) ≤ 1e-10

# With duplicate id parameters, wherse the parameters have the same value and therefore
# the model can be imported
path_SBML = joinpath(@__DIR__, "Models", "brusselator.xml")
prn, cb = load_SBML(path_SBML; inline_kineticlaw_parameters = false)
@test string.(parameters(prn.rn)) == ["B", "A", "k1", "C"]
@test string(prn.p[3]) == "k1 => 1.0"

# With duplicate id parameters, wherse the parameters have the same value and therefore
# the model can be imported
path_SBML = joinpath(@__DIR__, "Models", "brusselator_kinetic_parameters_fail.xml")
@test_throws SBMLImporter.SBMLSupport begin
    prn, cb = load_SBML(path_SBML; inline_kineticlaw_parameters = false)
end
