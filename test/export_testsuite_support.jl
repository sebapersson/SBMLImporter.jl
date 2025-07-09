using CSV, DataFrames, Downloads

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "testsuite_support.jl"))

df_export = DataFrame()

# Semantic test-suite
_testcases_not_supported = get_semantic_not_pass()
testcases_not_supported = reduce(vcat, values(_testcases_not_supported))
for i in 1:1821
    i % 50 == 0 && @info "i = $i"
    test_case = repeat("0", 5 - length(string(i))) * string(i)
    test_case in testcases_not_supported && continue
    tags = get_test_tags(test_case)
    df_tmp = DataFrame(testCase = test_case,
                       componentTags = prod(tags[:componentTags] .* ";")[1:end-1],
                       testTags = prod(tags[:testTags] .* ";")[1:end-1],
                       testSuite = "semantic")
    df_export = vcat(df_export, df_tmp)
end

# Stochastic test-suite
_testcases_not_supported = get_stochastic_not_pass()
testcases_not_supported = reduce(vcat, values(_testcases_not_supported))
for i in 1:100
    i % 10 == 0 && @info "i = $i"
    test_case = repeat("0", 5 - length(string(i))) * string(i)
    test_case in testcases_not_supported && continue
    tags = get_test_tags(test_case)
    df_tmp = DataFrame(testCase = test_case,
                       componentTags = prod(tags[:componentTags] .* ";")[1:end-1],
                       testTags = prod(tags[:testTags] .* ";")[1:end-1],
                       testSuite = "stochastic")
    df_export = vcat(df_export, df_tmp)
end

CSV.write(joinpath(@__DIR__, "..", "sbml_testsuite_support.csv"))
