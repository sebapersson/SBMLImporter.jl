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

sbml_switch(cond) = """
<sbml level="3" version="2">
  <model id="switch">
    <listOfParameters><parameter id="x" value="0" constant="false"/></listOfParameters>
    <listOfRules><rateRule variable="x"><math>
      <piecewise><piece><cn>1</cn>$cond</piece><otherwise><cn>0</cn></otherwise></piecewise>
    </math></rateRule></listOfRules>
  </model>
</sbml>
"""
TIME = "<ci>time</ci>"
conditions = [
    "<apply><gt/>$TIME<cn>5</cn></apply>",                                     # t > 5
    "<apply><lt/>$TIME<cn>5</cn></apply>",                                     # t < 5
    "<apply><lt/><cn>5</cn>$TIME</apply>",                                     # 5 < t
    "<apply><gt/><apply><minus/>$TIME</apply><cn>-5</cn></apply>",             # -t > -5
    "<apply><gt/><apply><minus/>$TIME<cn>5</cn></apply><cn>0</cn></apply>",    # t - 5 > 0
    "<apply><lt/><apply><minus/>$TIME<cn>5</cn></apply><cn>0</cn></apply>",    # t - 5 < 0
    "<apply><gt/><apply><minus/><cn>5</cn>$TIME</apply><cn>0</cn></apply>",    # 5 - t > 0
]
for cond in conditions
    rn, cb = load_SBML(sbml_switch(cond); model_as_string = true)
    prob = ODEProblem(rn, get_u0_map(rn), (0.0, 10.0), get_parameter_map(rn))
    sol = solve(prob, Rodas5P(); callback = cb, tstops = [5.0])
    @test sol[rn.x][end] ≈ 5.0 atol = 1.0e-4
end
