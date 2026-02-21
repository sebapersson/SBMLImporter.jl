"""
    getcompartment(s)

For an SBML specie `s`, get its associated SBML `compartment`.

If `s` does not have a compartment, or is not a specie, `nothing` is returned.

## Example
```julia
# Get compartment for first model specie
using SBMLImporter, Catalyst
rn, cb = load_SBML(path_SBML)
sbml_species = species(rn)
c = getcompartment(sbml_species[1])
```
"""
function getcompartment(s)::Union{Nothing, String}
    @assert !(s isa AbstractVector) "Input s = $s must be a Catalyst specie (not a Vector)"
    md = ModelingToolkitBase.metadata(s)
    !haskey(md, CompartmentSBML) && return nothing
    return ModelingToolkitBase.getmetadata(s, CompartmentSBML)
end
