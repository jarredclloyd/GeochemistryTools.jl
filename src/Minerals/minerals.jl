export element_to_oxide, cation_to_oxide
# abstract type Mineral end

# struct Mica <: Mineral

# end

"""
    element_to_oxide(value; kwargs...)

Converts a value representing an amount of a given element into its corresponding oxide.

See also: `cation_to_oxide` for converting complex cations to oxides.

# Arguments:
- `value::Real`: The amount of the element to be converted, set to 1 if you just want the oxide conversion factor.

# Keyword Arguments:
- `element::Union{AbstractString, AbstractChar, Symbol, Integer}`: A string (`"Si"`), character (`'H'`), `Symbol` (`:Na`), or integer (`1`) indicating the element being converted.
- `element_multiplicity::Integer`: The stoichiometric number of the element present the oxide.
    + For example, for Al2O3, set this to `2`
    + default is `1`.
- `oxide::AbstractString`: The formula of the oxide to which the element should be converted.
    + For example, Al2O3.
- `units::AbstractString`: The units of the output value ("wt%" or "ppm").
    + Default is "wt%".

# Examples
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

"""
    cation_to_oxide(value; kwargs...)

Converts a value representing an amount of a given cation into its corresponding oxide.

See also: `element_to_oxide` for simple element conversions.

# Arguments:
- `value::Real`: The amount of the cation to be converted, set to 1 if you just want the oxide conversion factor.

# Keyword Arguments:
- `cation::Union{AbstractString}`: A string (`"NH4"`) indicating the cation being converted.
- `cation_multiplicity::Integer`: The stoichiometric number of the cation present the oxide.
    + For example, for "(NH4)2O", set this to `2`
    + default is `1`.
- `oxide::AbstractString`: The formula of the oxide to which the cation should be converted.
    + For example "(NH4)2O"
- `units::AbstractString`: The units of the output value ("wt%" or "ppm").
    + Default is "wt%".

# Examples
```julia-repl
julia> cation_to_oxide(1; cation=:"NH4", cation_multiplicity = 2, oxide = "(NH4)2O")
1.44345584566772
```
"""
function cation_to_oxide(value::Real; cation::Union{AbstractString}, cation_multiplicity::Integer=1, oxide::AbstractString, units::AbstractString="wt%")
    if units == "wt%"
        return value / (molecular_mass(cation; verbose=false) * cation_multiplicity) * molecular_mass(oxide; verbose=false)
    elseif units == "ppm"
        return (value * 1e-4) / (molecular_mass(cation; verbose=false) * cation_multiplicity) * molecular_mass(oxide; verbose=false)
    else
        throw(ArgumentError("""Invalid units. Please choose `"wt%"`  or `"ppm"`."""))
    end
end
