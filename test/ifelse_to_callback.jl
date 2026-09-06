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

# Time can appear in many forms in a piecewise condition, and each form must be rewritten
# to a callback that switches at the correct time, and that activates the correct piecewise
# branch. Historically only conditions on the form t - c {>, ≥, <, ≤} 0 were correctly
# handled, see https://github.com/sebapersson/SBMLImporter.jl/issues/168
function solve_sbml(
        path::String; ifelse_to_callback::Bool = true, model_as_string::Bool = false,
        tspan = (0.0, 10.0), pset = Dict{Symbol, Float64}(), kwargs...
    )
    rn, cb = load_SBML(
        path; ifelse_to_callback = ifelse_to_callback, model_as_string = model_as_string
    )
    sys = mtkcompile(ode_model(rn))
    vmap = merge(Dict(get_u0_map(rn)), Dict(get_parameter_map(rn)))
    for (id, value) in pset
        vmap[getproperty(rn, id)] = value
    end
    oprob = ODEProblem(sys, vmap, tspan)
    sol = solve(
        oprob, Rodas5P(); callback = cb, abstol = 1.0e-10, reltol = 1.0e-10, kwargs...
    )
    return sys, sol
end

const TIME_MATHML = "<csymbol encoding=\"text\" definitionURL=\"http://www.sbml.org/sbml/symbols/time\"> time </csymbol>"
_cn(value) = "<cn> $(value) </cn>"
_apply(fn::String, args...) = "<apply><$(fn)/>" * prod(args) * "</apply>"

# dx/dt = 1 if the condition holds and dx/dt = 2 otherwise. As x(0) = 0 and every
# condition below switches value at t = 2, x is given by the analytical solutions below
function _model_piecewise_time(condition::String)::String
    return """<?xml version="1.0" encoding="UTF-8"?>
    <sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
      <model id="piecewise_time">
        <listOfParameters>
          <parameter id="x" value="0" constant="false"/>
        </listOfParameters>
        <listOfRules>
          <rateRule variable="x">
            <math xmlns="http://www.w3.org/1998/Math/MathML">
              <piecewise>
                <piece>$(_cn(1))$(condition)</piece>
                <otherwise>$(_cn(2))</otherwise>
              </piecewise>
            </math>
          </rateRule>
        </listOfRules>
      </model>
    </sbml>
    """
end
_x_holds_before(t) = t ≤ 2.0 ? t : 2.0 + 2.0 * (t - 2.0)
_x_holds_after(t) = t ≤ 2.0 ? 2.0 * t : 4.0 + (t - 2.0)

const T = TIME_MATHML
for (fn, holds) in [("gt", :after), ("geq", :after), ("lt", :before), ("leq", :before)]
    # holds is when the condition t fn 2 holds. If the arguments are reversed, or if the
    # expression with time decreases with time, the opposite branch holds with time
    x_time_lhs = holds == :after ? _x_holds_after : _x_holds_before
    x_time_rhs = holds == :after ? _x_holds_before : _x_holds_after
    conditions = [
        # Time compared directly against a constant
        ("t $fn 2", _apply(fn, T, _cn(2)), x_time_lhs),
        ("2 $fn t", _apply(fn, _cn(2), T), x_time_rhs),
        # Time in a shifted expression
        ("t - 2 $fn 0", _apply(fn, _apply("minus", T, _cn(2)), _cn(0)), x_time_lhs),
        ("0 $fn t - 2", _apply(fn, _cn(0), _apply("minus", T, _cn(2))), x_time_rhs),
        ("2 - t $fn 0", _apply(fn, _apply("minus", _cn(2), T), _cn(0)), x_time_rhs),
        ("0 $fn 2 - t", _apply(fn, _cn(0), _apply("minus", _cn(2), T)), x_time_lhs),
        # Negated time
        ("-t $fn -2", _apply(fn, _apply("minus", T), _cn(-2)), x_time_rhs),
        # Time scaled by a constant
        ("2t $fn 4", _apply(fn, _apply("times", _cn(2), T), _cn(4)), x_time_lhs),
        ("t / 2 $fn 1", _apply(fn, _apply("divide", T, _cn(2)), _cn(1)), x_time_lhs),
        ("-2t $fn -4", _apply(fn, _apply("times", _cn(-2), T), _cn(-4)), x_time_rhs),
        # Time scaled by a constant in a shifted expression
        (
            "2t - 4 $fn 0",
            _apply(fn, _apply("minus", _apply("times", _cn(2), T), _cn(4)), _cn(0)),
            x_time_lhs,
        ),
    ]
    for (name, condition, x_expected) in conditions
        @testset "Piecewise with $name in condition" begin
            _, sol = solve_sbml(
                _model_piecewise_time(condition); model_as_string = true,
                saveat = 0.0:0.5:10.0
            )
            @test all(.≈(sol[:x], x_expected.(sol.t), atol = 1.0e-6))
        end
    end
end

# The Weber model has piecewise expressions where time is compared directly against a
# constant (t < 0), and where time appears in a shifted expression (t - PdBu_time < 0)
@testset "Piecewise with time in condition Weber model" begin
    path_SBML = joinpath(@__DIR__, "Models", "model_Weber_BMC2015.xml")
    # Default parameter values switch off every piecewise branch, hence values that
    # trigger the piecewise expressions within the time-span are used
    pset = Dict(
        :Ect_Expr_PI4K3beta_flag => 1.0, :Ect_Expr_CERT_flag => 1.0, :PdBu_dose => 2.0,
        :PdBu_time => 30.0, :kb_NB142_70_dose => 3.0, :kb_NB142_70_time => 60.0
    )
    kwargs = (
        pset = pset, tspan = (0.0, 100.0), saveat = 0.0:5.0:100.0,
        tstops = [30.0, 60.0],
    )
    sys1, sol1 = solve_sbml(path_SBML; kwargs...)
    _, sol2 = solve_sbml(path_SBML; ifelse_to_callback = false, kwargs...)
    @test sol1.retcode == ReturnCode.Success
    # u3-u6 are the piecewise assignment rule variables, and they all appear in reaction
    # kinetic laws, hence any error in the rewritten piecewise propagates to the species
    for name in unknowns(sys1)
        @test all(.≈(sol1[name], sol2[name], rtol = 1.0e-6, atol = 1.0e-8))
    end
end
