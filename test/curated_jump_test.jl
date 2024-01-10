# Compares a model (with non-unitary stoichiometries) for which we know the catalyst reaction system it should be produced.
# Checks that simulations are correct and generates mass action jumps. 

# Fetch packages.
using Catalyst, SBMLImporter, Test

# Creates models.
rn_sbml, u0_sbml, ps_sbml = SBML_to_ReactionSystem("Models/brusselator.xml")

rn_catalyst = @reaction_network $(rn_sbml.name) begin
    @species Y(t) X(t) # SBMLImporter has flipped order of species and parameters.
    @parameters B A
    A*compartment, ∅ --> X
    2*compartment, 2X + Y --> 3X
    B*compartment, X --> Y
    1*compartment, X --> ∅
end
u0_catalyst = [:X => 2.0, :Y => 10.0]
ps_catalyst = [:B => 4.0, :A => 1.0, :compartment => 1.0]

# Tests model properties.
@test issetequal(species(rn_sbml), species(rn_catalyst))
@test issetequal(parameters(rn_sbml), parameters(rn_catalyst))
@test isequal(reactions(rn_sbml)[1], reactions(rn_catalyst)[1])
@test isequal(reactions(rn_sbml)[2], reactions(rn_catalyst)[2]) # Fails for weird Symbolics reason were (when loaded) 2*compartment != compartment / (1//2). Makes no sense. Evaluates to teh same though.
@test isequal(reactions(rn_sbml)[3], reactions(rn_catalyst)[3])
@test isequal(reactions(rn_sbml)[4], reactions(rn_catalyst)[4])

# Makes and tests jump simulations.
dprob_sbml = DiscreteProblem(rn_sbml, u0_sbml, (0.0,100.0), ps_sbml)
jprob_sbml = JumpProblem(rn_sbml, dprob_sbml, Direct())
sol_sbml = solve(jprob_sbml, SSAStepper(); seed=1234)

dprob_catalyst = DiscreteProblem(rn_catalyst, u0_catalyst, (0.0,100.0), ps_catalyst)
jprob_catalyst = JumpProblem(rn_catalyst, dprob_catalyst, Direct())
sol_catalyst = solve(jprob_catalyst, SSAStepper(); seed=1234)

@test sol_sbml == sol_catalyst

# Tests that both generates mass actions jumps.
@test jprob_sbml.massaction_jump.scaled_rates == jprob_catalyst.massaction_jump.scaled_rates
@test jprob_sbml.massaction_jump.reactant_stoch == jprob_catalyst.massaction_jump.reactant_stoch
@test jprob_sbml.massaction_jump.net_stoch == jprob_catalyst.massaction_jump.net_stoch