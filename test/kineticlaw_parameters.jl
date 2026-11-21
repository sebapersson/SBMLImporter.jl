using Catalyst, OrdinaryDiffEqRosenbrock, SBMLImporter, Test

#=
    Test that if kinetic law parameters are inlined correctly if the user selects that
    option. Further, tests different duplicate id cases.
=#

# First case with no duplicate ids. Parameter order is not guaranteed (depends on Dict
# iteration order, which is not stable across Julia versions), so parameters are checked
# by name, not position.
path_SBML = joinpath(@__DIR__, "Models", "brusselator_kinetic_parameters.xml")
rn1, cb1 = load_SBML(path_SBML; inline_kineticlaw_parameters = :promote)
u01, ps1 = get_u0_map(rn1), get_parameter_map(rn1)
@test issetequal(string.(parameters(rn1)), ["B", "A", "k1_r2", "k1_r4", "compartment"])
ps1_by_name = Dict(string(p) => v for (p, v) in ps1)
@test ps1_by_name["k1_r2"] == 1.0
@test ps1_by_name["k1_r4"] == 1.0

rn2, cb2 = load_SBML(path_SBML; inline_kineticlaw_parameters = :inline)
u02, ps2 = get_u0_map(rn2), get_parameter_map(rn2)
@test issetequal(string.(parameters(rn2)), ["B", "A", "compartment"])
oprob1 = ODEProblem(rn1, u01, (0.0, 10.0), ps1)
sol1 = solve(oprob1, Rodas5P())
oprob2 = ODEProblem(rn2, u02, (0.0, 10.0), ps2)
sol2 = solve(oprob2, Rodas5P())
@test sum(reduce(vcat, sol1.u .- sol2.u) .^ 2) ≤ 1.0e-10

# With duplicate id parameters, where the parameters have the same value and therefore
# the model can be imported
path_SBML = joinpath(@__DIR__, "Models", "brusselator.xml")
rn, cb = load_SBML(path_SBML; inline_kineticlaw_parameters = :promote)
@test issetequal(string.(parameters(rn)), ["B", "A", "k1", "C"])
ps_by_name = Dict(string(p) => v for (p, v) in get_parameter_map(rn))
@test ps_by_name["k1"] == 1.0

# With duplicate id parameters, where the parameters have the same value and therefore
# the model cannot be imported
path_SBML = joinpath(@__DIR__, "Models", "brusselator_kinetic_parameters_fail.xml")
@test_throws SBMLImporter.SBMLSupport begin
    _, _ = load_SBML(path_SBML; inline_kineticlaw_parameters = :promote)
end

# Repeated kinetic law ids with different values (issue #178). Not possible to promote
# with the SBML ids, but possible when promoting with qualified (reaction id prefixed) ids
path_SBML = joinpath(@__DIR__, "Models", "repeated_kinetic_parameters.xml")
@test_throws SBMLImporter.SBMLSupport begin
    _, _ = load_SBML(path_SBML; inline_kineticlaw_parameters = :promote)
end

rn_qualified, cb_qualified = load_SBML(
    path_SBML; inline_kineticlaw_parameters = :promote_qualified
)
ps_qualified = get_parameter_map(rn_qualified)
@test issetequal(string.(parameters(rn_qualified)), ["r1_k1", "r2_k1", "cell"])
ps_qualified_by_name = Dict(string(p) => v for (p, v) in ps_qualified)
@test ps_qualified_by_name["r1_k1"] == 0.9
@test ps_qualified_by_name["r2_k1"] == 3.6e-5

# Same dynamics as when the kinetic law parameters are inlined
rn_inline, cb_inline = load_SBML(path_SBML; inline_kineticlaw_parameters = :inline)
@test issetequal(string.(parameters(rn_inline)), ["cell"])
oprob_qualified = ODEProblem(
    rn_qualified, get_u0_map(rn_qualified), (0.0, 10.0), ps_qualified
)
oprob_inline = ODEProblem(
    rn_inline, get_u0_map(rn_inline), (0.0, 10.0), get_parameter_map(rn_inline)
)
sol_qualified = solve(
    oprob_qualified, Rodas5P(), abstol = 1.0e-10, reltol = 1.0e-10, saveat = 0:1:10
)
sol_inline = solve(
    oprob_inline, Rodas5P(), abstol = 1.0e-10, reltol = 1.0e-10, saveat = 0:1:10
)
@test sum(reduce(vcat, sol_qualified.u .- sol_inline.u) .^ 2) ≤ 1.0e-10

# Invalid option
@test_throws ArgumentError begin
    _, _ = load_SBML(path_SBML; inline_kineticlaw_parameters = :not_an_option)
end

# Providing a Bool is deprecated, but must give the same result as the corresponding Symbol
path_SBML = joinpath(@__DIR__, "Models", "brusselator_kinetic_parameters.xml")
rn_true = @test_deprecated r"inline_kineticlaw_parameters" first(
    load_SBML(path_SBML; inline_kineticlaw_parameters = true)
)
rn_false = @test_deprecated r"inline_kineticlaw_parameters" first(
    load_SBML(path_SBML; inline_kineticlaw_parameters = false)
)
@test issetequal(string.(parameters(rn_true)), string.(parameters(rn2)))
@test issetequal(string.(parameters(rn_false)), string.(parameters(rn1)))
