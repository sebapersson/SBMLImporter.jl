# Compares a model (with non-unitary stoichiometries) for which we know the catalyst
# reaction system it should be produced. Checks that simulations are correct and generates
# mass action jumps.

using Catalyst, JumpProcesses, SBMLImporter, OrdinaryDiffEqRosenbrock, Test

# Creates models.
rn_sbml, cb = load_SBML(joinpath(@__DIR__, "Models", "brusselator.xml"), massaction = true)
u0_sbml = get_u0_map(rn_sbml)
ps_sbml = get_parameter_map(rn_sbml)

rn_catalyst = @reaction_network $(nameof(rn_sbml)) begin
    @species Y(t) X(t) # SBMLImporter has flipped order of species and parameters.
    @parameters B A
    A * C, ∅ --> X
    2 * C, 2X + Y --> 3X
    B * C, X --> Y
    1 * C, X --> ∅
end
rn_catalyst = Catalyst.complete(rn_catalyst)
u0_catalyst = [:X => 2.0, :Y => 10.0]
ps_catalyst = [:B => 4.0, :A => 1.0, :C => 1.0]

# Tests model properties.
@test issetequal(species(rn_sbml), species(rn_catalyst))
@test issetequal(parameters(rn_sbml), parameters(rn_catalyst))
@test isequal(reactions(rn_sbml)[1], reactions(rn_catalyst)[1])
@test isequal(reactions(rn_sbml)[2], reactions(rn_catalyst)[4])
# @test isequal(reactions(rn_sbml)[3], reactions(rn_catalyst)[2]) # Mathematically equal
@test isequal(reactions(rn_sbml)[4], reactions(rn_catalyst)[3])

# Makes and tests jump simulations. Note that tests need to account for reactions not
# appearing in the same order
jprob_sbml = JumpProblem(rn_sbml, u0_sbml, (0.0, 100.0), ps_sbml)
jprob_catalyst = JumpProblem(rn_catalyst, u0_catalyst, (0.0, 100.0), ps_catalyst)
imap = [1, 4, 2, 3]
@test jprob_sbml.massaction_jump.scaled_rates == jprob_catalyst.massaction_jump.scaled_rates[imap]
@test jprob_sbml.massaction_jump.reactant_stoch == jprob_catalyst.massaction_jump.reactant_stoch[imap]
@test jprob_sbml.massaction_jump.net_stoch == jprob_catalyst.massaction_jump.net_stoch[imap]

# Check consistent simulations when model is imported with and without rewriting to
# Catalyst mass-action format
path = joinpath(@__DIR__, "Models", "brusselator.xml")
rn1, cb1 = load_SBML(path; massaction = true)
rn2, cb2 = load_SBML(path; massaction = false)
u01, ps1 = get_u0_map(rn1), get_parameter_map(rn1)
u02, ps2 = get_u0_map(rn2), get_parameter_map(rn2)
oprob1 = ODEProblem(rn1, u01, (0.0, 10.0), ps1)
oprob2 = ODEProblem(rn2, u02, (0.0, 10.0), ps2)
sol1 = solve(oprob1, Rodas5())
sol2 = solve(oprob2, Rodas5())
@test sol1 == sol2
