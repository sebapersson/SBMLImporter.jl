using SBMLImporter, Catalyst, Test

path_SBML = joinpath(@__DIR__, "Models", "brusselator.xml")
prn, cb = load_SBML(path_SBML)
sbml_reactions = reactions(prn.rn)
@test all(getmetadata.(sbml_reactions, :id) .== ["r1_id", "r4_id", "r2_id", "r3_id"])
@test all(getmetadata.(sbml_reactions, :name) .== ["r1_name", "r4_name", "r2_name", "r3_name"])
# Need to check that meta-data carries through if massaction is enabled
prn, cb = load_SBML(path_SBML; massaction=true)
sbml_reactions = reactions(prn.rn)
@test all(getmetadata.(sbml_reactions, :id) .== ["r1_id", "r4_id", "r2_id", "r3_id"])
@test all(getmetadata.(sbml_reactions, :name) .== ["r1_name", "r4_name", "r2_name", "r3_name"])
