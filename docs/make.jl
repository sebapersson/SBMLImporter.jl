using SBMLImporter
using Documenter

DocMeta.setdocmeta!(SBMLImporter, :DocTestSetup, :(using SBMLImporter); recursive = true)

format = Documenter.HTML(;prettyurls = get(ENV, "CI", "false") == "true",
                         assets = String["assets/custom_theme.css"],
                         repolink = "https://github.com/sebapersson/SBMLImporter.jl",
                         edit_link = "main")

makedocs(; modules = [SBMLImporter],
         repo = "https://github.com/sebapersson/SBMLImporter.jl/blob/{commit}{path}#{line}",
         checkdocs = :exports,
         warnonly = false,
         format = format,
         sitename = "SBMLImporter.jl",
         pages = [
             "Home" => "index.md",
             "Tutorial" => "tutorial.md",
             "API" => "API.md",
             "Supported SBML features" => "support.md",
             "Other SBML related Julia packages" => "differences.md",
             "FAQs" => "FAQ.md",
         ],)

deploydocs(; repo = "github.com/sebapersson/SBMLImporter.jl.git",
           devbranch = "main",)
