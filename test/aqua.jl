using Aqua
using SBMLImporter

@testset "Aqua" begin
    Aqua.test_ambiguities(SBMLImporter, recursive = false)
    Aqua.test_undefined_exports(SBMLImporter)
    Aqua.test_unbound_args(SBMLImporter)
    Aqua.test_stale_deps(SBMLImporter)
    Aqua.test_deps_compat(SBMLImporter)
    Aqua.find_persistent_tasks_deps(SBMLImporter)
    Aqua.test_piracies(SBMLImporter)
    Aqua.test_project_extras(SBMLImporter)
end
