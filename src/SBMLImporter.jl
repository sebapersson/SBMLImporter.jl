module SBMLImporter

using Catalyst: Catalyst, setmetadata, parameters, unknowns, @unpack, get_u0_map,
    get_parameter_map
using ComponentArrays: ComponentArray
using DiffEqBase: CallbackSet, DiscreteCallback, ContinuousCallback
using JumpProcesses: reset_aggregated_jumps!
using PrecompileTools: @setup_workload, @compile_workload
import ModelingToolkitBase
using RuntimeGeneratedFunctions: RuntimeGeneratedFunctions, @RuntimeGeneratedFunction
using SBML: SBML, readSBMLFromString
import SciMLBase
using Setfield: @set
import SpecialFunctions
import Symbolics

RuntimeGeneratedFunctions.init(@__MODULE__)

include("structs.jl")

const SBMLMathVariables = Union{
    SBML.MathIdent, SBML.MathVal, SBML.MathTime, SBML.MathConst, SBML.MathAvogadro,
}
const SBMLRule = Union{SBML.AssignmentRule, SBML.RateRule, SBML.AlgebraicRule}
const FORBIDDEN_IDS = ["true", "false", "time", "pi", "Inf", "NaN", "Differential"]
const VariableSBML = Union{SpecieSBML, ParameterSBML, CompartmentSBML}

include("callbacks.jl")
include("common.jl")
include("compartments.jl")
include("events.jl")
include("functions.jl")
include("initial_assignments.jl")
include("load.jl")
include("math.jl")
include("odeproblem.jl")
include("parameters.jl")
include("parse.jl")
include("piecewise.jl")
include("reactions.jl")
include("replace_idents.jl")
include("rules.jl")
include("sbml_functions.jl")
include("support.jl")
include("species.jl")
include("system.jl")
include("templates.jl")
include("util.jl")

@setup_workload begin
    dirmodels = joinpath(@__DIR__, "..", "test", "Models")
    # Model without events
    path_no_events = joinpath(dirmodels, "model_Boehm_JProteomeRes2014.xml")
    # Model with events
    path_events = joinpath(dirmodels, "model_Brannmark_JBC2010.xml")
    @compile_workload begin
        rn, cb = load_SBML(path_no_events)
        rn, cb = load_SBML(path_events)
    end
end

export load_SBML, getcompartment, get_u0_map, get_parameter_map

end
