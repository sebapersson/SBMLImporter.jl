module SBMLImporterForwardDiffExt

import ForwardDiff
import SBMLImporter

"""
    dual_to_float(x::ForwardDiff.Dual)::AbstractFloat

Via recursion convert a Dual to a Float.
"""
function SBMLImporter._to_float(x::ForwardDiff.Dual)::AbstractFloat
    return SBMLImporter._to_float(x.value)
end

end
