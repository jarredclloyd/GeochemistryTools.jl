export element_to_oxide
# abstract type Mineral end

# struct Mica <: Mineral

# end

"""
    element_to_oxide(value; kwargs...)

Converts a value representing an amount of a given element into its corresponding oxide.

# Arguments:
- `value::Real`: The amount of the element to be converted, set to 1 if you just want the oxide conversion factor.

# Keyword Arguments:
- `element::Union{AbstractString, AbstractChar}`: A string ("Si") or character ('H') indicating the element being converted.
- `element_multiplicity::Integer`: The stoichiometric number of the element present the oxide.
    + For example, for Al2O3, set this to `2`
    + default is `1`.
- `oxide::AbstractString`: The formula of the oxide to which the element should be converted.
- `units::AbstractString`: The units of the output value ("wt%" or "ppm").
    + Default is "wt%".



```julia-repl
julia> element_to_oxide(1; element="Si",element_multiplicity=1,oxide="SiO2")
2.1393352441651383
```
"""
function element_to_oxide(value::Real; element::Union{AbstractString, AbstractChar, Symbol, Integer}, element_multiplicity::Integer=1, oxide::AbstractString, units::AbstractString="wt%")
    if units == "wt%"
        return value / (get_atomicmass(element) * element_multiplicity) * molecular_mass(oxide; verbose=false)
    elseif units == "ppm"
        return (value * 1e-4) / (get_atomicmass(element) * element_multiplicity) * molecular_mass(oxide; verbose=false)
    else
        throw(ArgumentError("""Invalid units. Please choose `"wt%"`  or `"ppm"`."""))
    end
end
