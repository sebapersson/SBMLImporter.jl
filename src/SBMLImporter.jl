module SBMLImporter

using Catalyst
using DiffEqBase
using DiffEqCallbacks
using JumpProcesses
using PrecompileTools
using ReactionNetworkImporters
using RuntimeGeneratedFunctions
using SBML
using Setfield
using SpecialFunctions

RuntimeGeneratedFunctions.init(@__MODULE__)

# Importing SBML models
include("Structs.jl")
include("Build_model.jl")
include("Callbacks.jl")
include("Check_support.jl")
include("SBML_to_ReactionSystem.jl")
include("Parameters.jl")
include("Species.jl")
include("Functions.jl")
include("Rules.jl")
include("Common.jl")
include("Piecewise.jl")
include("Events.jl")
include("Math.jl")
include("Initial_assignments.jl")
include("Reactions.jl")

@setup_workload begin
    # Model without events 
    path_SBML = joinpath(@__DIR__, "..", "test", "Models",
                         "model_Boehm_JProteomeRes2014.xml")
    parsed_rn, cb = load_SBML(path_SBML)
    # Model with events 
    path_SBML = joinpath(@__DIR__, "..", "test", "Models", "model_Brannmark_JBC2010.xml")
    parsed_rn, cb = load_SBML(path_SBML)
end

export load_SBML

end
