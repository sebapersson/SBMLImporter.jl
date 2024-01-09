#=
    Test that import time for a larger model of 1265 species and 1000+ reactions
=#

using SBMLImporter
using Downloads
using Test


#sbml_url = "https://raw.githubusercontent.com/SciML/Catalyst_PLOS_COMPBIO_2023/master/Benchmarks/Data/BCR_no_obs.xml"
sbml_url = "https://www.ebi.ac.uk/biomodels/model/download/MODEL1112100000.2?filename=MODEL1112100000_url.xml"
sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))

b1 = @elapsed rn, specie_map, parameter_map, cb = SBML_to_ReactionSystem(sbml_string; model_as_string=true, ret_all=true)

# Usually takes around 4s locally, but better to brace for GitHub CI
@test b1 ≤ 20
