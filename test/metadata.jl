using Catalyst, SBMLImporter, Test

path_SBML = joinpath(@__DIR__, "Models", "brusselator.xml")
rn, cb = load_SBML(path_SBML)
sbml_reactions = reactions(rn)
@test all(getmetadata.(sbml_reactions, :id) .== ["r1_id", "r4_id", "r2_id", "r3_id"])
@test all(getmetadata.(sbml_reactions, :name) .== ["r1_name", "r4_name", "r2_name", "r3_name"])

# Need to check that meta-data carries through if mass-action is enabled
rn, cb = load_SBML(path_SBML; massaction = true)
sbml_reactions = reactions(rn)
@test all(getmetadata.(sbml_reactions, :id) .== ["r1_id", "r4_id", "r2_id", "r3_id"])
@test all(getmetadata.(sbml_reactions, :name) .== ["r1_name", "r4_name", "r2_name", "r3_name"])
sbml_species = species(rn)
@test getcompartment.(sbml_species) == ["C", "C"]

# Test correct compartments for two compartment model
path_SBML = joinpath(@__DIR__, "Models", "model_Boehm_JProteomeRes2014.xml")
rn, cb = load_SBML(path_SBML)
sbml_species = species(rn)
@test getcompartment.(sbml_species) == ["cyt", "cyt", "nuc", "nuc", "cyt", "cyt", "nuc", "cyt"]
