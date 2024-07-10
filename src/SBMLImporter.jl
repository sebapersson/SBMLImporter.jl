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

include("Structs.jl")

const SBMLMathVariables = Union{SBML.MathIdent, SBML.MathVal, SBML.MathTime, SBML.MathConst, SBML.MathAvogadro}
const SBMLRule = Union{SBML.AssignmentRule, SBML.RateRule, SBML.AlgebraicRule}
const FORBIDDEN_IDS = ["true", "false", "time", "pi", "Inf", "NaN"]
const VariableSBML = Union{SpecieSBML, ParameterSBML, CompartmentSBML}

include("Build_model.jl")
include("Callbacks.jl")
include("compartments.jl")
include("Check_support.jl")
include("SBML_to_ReactionSystem.jl")
include("parameters.jl")
include("Species.jl")
include("Functions.jl")
include("Rules.jl")
include("Common.jl")
include("Piecewise.jl")
include("Events.jl")
include("sbml_functions.jl")
include("math.jl")
include("replace_idents.jl")
include("templates.jl")
include("system.jl")
include("Initial_assignments.jl")
include("Reactions.jl")

#=
@setup_workload begin
    dirmodels = joinpath(@__DIR__, "..", "test", "Models")
    # Model without events
    path_SBML = joinpath(dirmodels, "model_Boehm_JProteomeRes2014.xml")
    parsed_rn, cb = load_SBML(path_SBML)
    # Model with events
    path_SBML = joinpath(dirmodels, "model_Brannmark_JBC2010.xml")
    parsed_rn, cb = load_SBML(path_SBML)
end
=#

export load_SBML

end
