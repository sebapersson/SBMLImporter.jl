function get_semantic_not_pass()
    not_pass = Dict{Symbol, Vector{String}}()
    # Delay differential equations not yet supported
    not_pass[:delay] = [
        "00937", "00938", "00939", "00940", "00941", "00942", "00943",
        "00981", "00982", "00983", "00984", "00985", "01400", "01401",
        "01403", "01404", "01410", "01411", "01412", "01413", "01414",
        "01415", "01416", "01417", "01418", "01419", "01454", "01480",
        "01481", "01518", "01522", "01523", "01524", "01534", "01535",
        "01536", "01537", "01538", "01539", "01406", "01407", "01409",
        "01410",
    ]
    # Events with delay not yet supporter
    not_pass[:event_delay] = [
        "00071", "00072", "00073", "00405", "00405", "00406", "00407", "00409", "00410",
        "00411", "00412", "00413", "00414", "00415", "00416", "00417", "00418", "00419",
        "00420", "00622", "00623", "00624", "00637", "00638", "00639", "00649", "00650",
        "00651", "00664", "00665", "00666", "00682", "00683", "00684", "00690", "00702",
        "00708", "00724", "00737", "00757", "00758", "00759", "00763", "00764", "00765",
        "00766", "00767", "00768", "00769", "00770", "00771", "00772", "00773", "00774",
        "00775", "00776", "00777", "00778", "00779", "00780", "00848", "00849", "00850",
        "00886", "00887", "00932", "00933", "00936", "00408", "00461", "00655", "00656",
        "00657", "00980", "01000", "01048", "01049", "01050", "01074", "01075", "01076",
        "01119", "01120", "01230", "01241", "01263", "01268", "01269", "01270", "01287",
        "01295", "01299", "01305", "01324", "01325", "01326", "01327", "01328", "01329",
        "01335", "01507", "01508", "01509", "01511", "01519", "01520", "01525", "01526",
        "01528", "01529", "01532", "01575", "01576", "01579", "01580", "01581", "01582",
        "01584", "01585", "01586", "01587", "01594", "01595", "01597", "01598", "01600",
        "01601", "01602", "01603", "01604", "01659", "01660", "01661", "01673", "01674",
        "01675", "01676", "01677", "01678", "01679", "01680", "01687", "01688", "01689",
        "01690", "01691", "01692", "01701", "01702", "01703", "01704", "01706", "01707",
        "01708", "01709", "01710", "01711", "01712", "01713", "01715", "01716", "01717",
        "01718", "01720", "01721", "01754", "01755", "01756", "01757", "01758", "01759",
        "01798", "00731", "01095", "01771", "01672", "00701", "01305",
    ]
    not_pass[:event_delay] = vcat(not_pass[:event_delay], ["004" * string(i) for i in 21:61])
    # Events with priority not yet supported
    not_pass[:event_priority] = [
        "00931", "00934", "00935", "00962", "00963", "00964", "00965", "00966", "00967",
        "00978", "00978", "01229", "01242", "01267", "01294", "01298", "01331", "01332",
        "01333", "01334", "01336", "01337", "01466", "01512", "01521", "01533", "01577",
        "01583", "01588", "01589", "01590", "01591", "01592", "01593", "01599", "01605",
        "01626", "01627", "01662", "01681", "01682", "01683", "01705", "01714", "01212",
        "01213", "01214", "01215", "00952", "00930", "00997", "01262", "01286", "01330",
        "01772",
    ]
    # Fast reactions not yet supported
    not_pass[:fast_reactions] = [
        "00870", "00871", "00872", "00873", "00874", "00875", "00986", "00987", "00988",
        "01051", "01052", "01053", "01396", "01397", "01398", "01399", "01544", "01545",
        "01546", "01547", "01548", "01549", "01550", "01551", "01558", "01559", "01560",
        "01565", "01567", "01568", "01569", "01570", "01571", "01572",
    ]
    # Uncommon mat expressions not supported
    not_pass[:imply] = ["01274", "01279", "01497"]
    not_pass[:three_arg_and_xor] = ["00198", "00200", "00201", "00276", "00279", "01282"]
    not_pass[:piecewise_algebraic_rules] = ["01503"]
    not_pass[:rem_div] = ["01272", "01273", "01277", "01278", "01495", "01496"]
    not_pass[:event_multiple_triggers] = ["01211", "01531"]
    not_pass[:multiple_arg_inequality] = [
        "01216", "01494", "01781", "01782", "01783", "01209", "01210", "00278",
    ]
    # Hierarchical models extension not yet supported
    not_pass[:hierarchical_models] = vcat(
        ["011" * string(i) for i in 24:83], ["0" * string(i) for i in 1344:1394],
        ["0" * string(i) for i in 1467:1477], ["01778"]
    )
    # FBA models will always be outside the package scope
    not_pass[:fba_models] = vcat(
        ["01" * string(i) for i in 186:196], ["01" * string(i) for i in 606:625],
        ["01627", "01628", "01629", "01630"]
    )
    # Certain variable names conflicting with Julia reserved variables
    not_pass[:bad_names] = [
        "01810", "01811", "01812", "01813", "01814", "01815", "01816", "01817", "01818",
        "01819", "01820", "01821",
    ]
    # Rare features which SBMLImporter cannot catch
    not_captured = String[]
    # Edge-cases where I think there is a bug in SBML.jl
    not_captured = vcat(not_captured, ["01318", "01319", "01320"])
    # Uncommon mathML with boundary effects, not sure I agree with test-suite
    not_captured = vcat(not_captured, ["00959", "01488"])
    # A form of event competition
    not_captured = vcat(not_captured, ["00953"])
    # If user wants to add a random species, it must either be as a species,
    # initialAssignment, assignmentRule or by event, not by just random adding it to
    # equations. Cannot capture this error.
    not_captured = vcat(not_captured, ["00974"])
    # Bug in SBML.jl (parameter rateOf)
    not_captured = vcat(not_captured, ["01321", "01322"])
    # Even trigger cannot depend on algebraic rule
    not_captured = vcat(not_captured, ["01578"])
    # Inf or NaN in initial assignments
    not_captured = vcat(not_captured, ["00950", "00951"])
    # Subset of empty reactions. TODO: Add fix for
    not_captured = vcat(not_captured, ["01304", "01305", "01306"])
    not_pass[:not_captured] = not_captured
    return not_pass
end

function get_stochastic_not_pass()
    not_pass = Dict{Symbol, Vector{String}}()

    not_pass[:assignment_rule] = ["00019"]
    not_pass[:continuous_callback] = ["00033"]
    not_pass[:distributions] = vcat(["000$i" for i in 40:99], ["00100"])
    return not_pass
end
