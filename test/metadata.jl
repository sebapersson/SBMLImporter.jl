using Catalyst, SBMLImporter, Test

# Reaction/species order is not guaranteed (depends on Dict iteration order, which is not
# stable across Julia versions)
expected_name_for_id = Dict(
    "r1_id" => "r1_name", "r4_id" => "r4_name", "r2_id" => "r2_name", "r3_id" => "r3_name"
)

path_SBML = joinpath(@__DIR__, "Models", "brusselator.xml")
rn, cb = load_SBML(path_SBML)
sbml_reactions = reactions(rn)
ids = getmetadata.(sbml_reactions, :id)
names = getmetadata.(sbml_reactions, :name)
@test issetequal(ids, keys(expected_name_for_id))
@test all(name == expected_name_for_id[id] for (id, name) in zip(ids, names))

# Need to check that meta-data carries through if mass-action is enabled
rn, cb = load_SBML(path_SBML; massaction = true)
sbml_reactions = reactions(rn)
ids = getmetadata.(sbml_reactions, :id)
names = getmetadata.(sbml_reactions, :name)
@test issetequal(ids, keys(expected_name_for_id))
@test all(name == expected_name_for_id[id] for (id, name) in zip(ids, names))
sbml_species = species(rn)
@test getcompartment.(sbml_species) == ["C", "C"]

# Test correct compartments for two compartment model
path_SBML = joinpath(@__DIR__, "Models", "model_Boehm_JProteomeRes2014.xml")
rn, cb = load_SBML(path_SBML)
sbml_species = species(rn)
expected_compartment = Dict(
    "STAT5A(t)" => "cyt", "pApA(t)" => "cyt", "nucpApB(t)" => "nuc",
    "nucpBpB(t)" => "nuc", "STAT5B(t)" => "cyt", "pApB(t)" => "cyt",
    "nucpApA(t)" => "nuc", "pBpB(t)" => "cyt"
)
@test issetequal(string.(sbml_species), keys(expected_compartment))
@test all(
    getcompartment(sp) == expected_compartment[string(sp)] for sp in sbml_species
)
