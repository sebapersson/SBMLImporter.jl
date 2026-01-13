using SBMLImporter
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(SBMLImporter, :DocTestSetup, :(using SBMLImporter); recursive = true)

format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    assets = String["assets/custom_theme.css"],
    repolink = "https://github.com/sebapersson/SBMLImporter.jl",
    edit_link = "main"
)

makedocs(;
    modules = [SBMLImporter],
    sitename = "SBMLImporter.jl",
    repo = Remotes.GitHub("sebapersson", "SBMLImporter.jl"),
    authors = "Sebastian Persson, and contributors",
    checkdocs = :exports,
    warnonly = false,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/sebapersson/SBMLImporter.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "API" => "API.md",
        "SBML support and ecosystem" =>
            Any[
            "Supported SBML features" => "support.md",
            "Other SBML related Julia packages" => "differences.md",
        ],
        "Contributing" => "contributing.md",
    ],
)

DocumenterVitepress.deploydocs(
    repo = "github.com/sebapersson/SBMLImporter.jl.git",
    target = "build", # this is where Vitepress stores its output
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
