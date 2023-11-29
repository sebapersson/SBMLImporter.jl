module SBMLImporter

using DiffEqCallbacks
using ModelingToolkit
using OrdinaryDiffEq
using PrecompileTools
using RuntimeGeneratedFunctions
using SBML
using SpecialFunctions

RuntimeGeneratedFunctions.init(@__MODULE__)

# Importing SBML models
include("Structs.jl")
include("Callbacks.jl")
include("Check_support.jl")
include("SBML_to_ModellingToolkit.jl")
include("SBML_parameters.jl")
include("SBML_species.jl")
include("SBML_functions.jl")
include("SBML_rules.jl")
include("Common.jl")
include("Piecewise.jl")
include("Events.jl")
include("SBML_math.jl")
include("Initial_assignments.jl")
include("Reactions.jl")


# Model with piecewise for PrecompileTools to make model reading faster
@setup_workload begin
    path_SBML = joinpath(@__DIR__, "..", "test", "Models", "model_Boehm_JProteomeRes2014.xml")
    ode_system, specie_map, parameter_map, cb, get_tstops, ifelse_t0 = SBML_to_ODESystem(path_SBML, ret_all=true, write_to_file=false)
    ode_problem = ODEProblem(ode_system, specie_map, (0.0, 10.0), parameter_map, jac=true)
    for _f! in ifelse_t0
        _f!(ode_problem.u0, ode_problem.p)
    end
    tstops = get_tstops(ode_problem.u0, ode_problem.p)
    sol = solve(ode_problem, Rodas5P(), abstol=1e-8, reltol=1e-8, tstops=tstops, callback=cb)
end

export SBML_to_ODESystem

end