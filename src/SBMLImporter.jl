module SBMLImporter

using Catalyst
using DiffEqBase
using DiffEqCallbacks
using JumpProcesses
using PrecompileTools
using ReactionNetworkImporters
using RuntimeGeneratedFunctions
using SBML
using SciMLBase
using Setfield
using SpecialFunctions

RuntimeGeneratedFunctions.init(@__MODULE__)

include("structs.jl")

const SBMLMathVariables = Union{SBML.MathIdent, SBML.MathVal, SBML.MathTime, SBML.MathConst, SBML.MathAvogadro}
const SBMLRule = Union{SBML.AssignmentRule, SBML.RateRule, SBML.AlgebraicRule}
const FORBIDDEN_IDS = ["true", "false", "time", "pi", "Inf", "NaN"]
const VariableSBML = Union{SpecieSBML, ParameterSBML, CompartmentSBML}

include("callbacks.jl")
include("common.jl")
include("compartments.jl")
include("events.jl")
include("functions.jl")
include("initial_assignments.jl")
include("load.jl")
include("math.jl")
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

@setup_workload begin
    dirmodels = joinpath(@__DIR__, "..", "test", "Models")
    # Model without events
    path_SBML = joinpath(dirmodels, "model_Boehm_JProteomeRes2014.xml")
    parsed_rn, cb = load_SBML(path_SBML)
    # Model with events
    path_SBML = joinpath(dirmodels, "model_Brannmark_JBC2010.xml")
    parsed_rn, cb = load_SBML(path_SBML)
end

export load_SBML

end
