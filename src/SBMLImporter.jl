module SBMLImporter

using DiffEqCallbacks
using ModelingToolkit
using OrdinaryDiffEq
using PrecompileTools
using RuntimeGeneratedFunctions
using SBML
using SpecialFunctions
using Catalyst
using Setfield
using JumpProcesses

RuntimeGeneratedFunctions.init(@__MODULE__)

# Importing SBML models
include("Structs.jl")
include("Build_model.jl")
include("Callbacks.jl")
include("Check_support.jl")
include("SBML_to_ODESystem.jl")
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

# Model with piecewise for PrecompileTools to make model reading faster
@setup_workload begin
    path_SBML = joinpath(@__DIR__, "..", "test", "Models", "model_Boehm_JProteomeRes2014.xml")
    ode_system, specie_map, parameter_map, cb = SBML_to_ODESystem(path_SBML, ret_all=true)
    ode_problem = ODEProblem(ode_system, specie_map, (0.0, 10.0), parameter_map, jac=true)
    sol = solve(ode_problem, Rodas5P(), abstol=1e-8, reltol=1e-8, callback=cb)
end

export SBML_to_ODESystem, SBML_to_ReactionSystem

end