#=
    Test that import time for a larger model of 1265 species and 1000+ reactions
=#

using SBMLImporter
using Downloads
using Test


sbml_url = "https://www.ebi.ac.uk/biomodels/model/download/MODEL1112100000.2?filename=MODEL1112100000_url.xml"
sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))

b1 = @elapsed parsed_rn, cb = load_SBML(sbml_string; model_as_string=true)

# Usually takes around 4s locally, but better to brace for GitHub CI
@test b1 ≤ 20
