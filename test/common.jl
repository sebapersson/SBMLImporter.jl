function get_sbml_urls(base_url::String)
    levels = ["2", "3"]
    sublevels = ["1", "2", "3", "4", "5"]
    sbml_urls = String[]
    for level in levels
        for sublevel in sublevels
            sbml_url = base_url * "-sbml-l" * level * "v" * sublevel * ".xml"
            try
                sbml = String(take!(Downloads.download(sbml_url, IOBuffer())))
                push!(sbml_urls, sbml_url)
            catch
            end
        end
    end
    return sbml_urls
end

# Function to get model-str
function get_model_str(test_case)
    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/semantic/$test_case/$test_case"
    sbml_urls = get_sbml_urls(base_url)
    sbml_url = sbml_urls[end]
    sbml_string = String(take!(Downloads.download(sbml_url, IOBuffer())))
    model_SBML = SBMLImporter.build_SBML_model(sbml_string, model_as_string = true,
                                               inline_assignment_rules = false)

    parsed_model_SBML = SBMLImporter._reactionsystem_from_SBML(model_SBML)
    model_str = SBMLImporter.reactionsystem_to_string(parsed_model_SBML, false, "",
                                                      model_SBML)

    return model_str, model_SBML
end
