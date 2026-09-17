# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

# GeochemistryToolsCore.jl/AccessorFunctions.jl
# Helper functions to access component of package dependencies.
export ionic_radius, atomic_mass

"""
    ionic_radius(element::Union{String,Symbol};
    coordination::Symbol, charge::Integer, [value_only::Bool = false])

Helper function to return the ionic radius (pm) of an element with a given charge state and coordination number.

Wraps around Mendeleev.elements[:Symbol], see [Mendeleev.jl](https://eben60.github.io/Mendeleev.jl/).

User must specify the coordination state as a roman numeral (string or symbol) and the charge state as a signed integer (e.g. -1, 1).

# Example

```julia-repl
julia> ionic_radius(:Nd; coordination = :VIII, charge = 3)
110.9 pm

julia> ionic_radius(\"Nd\"; coordination = :VIII, charge = 3, value_only = true)
110.9
```
"""
function ionic_radius(
    element::Union{String,Symbol,Int};
    coordination::Union{String,Symbol},
    charge::Int,
    value_only::Bool = false
)
    if isa(element, AbstractString)
        element = Symbol(element)
    end
    if isa(coordination, AbstractString)
        coordination = Symbol(coordination)
    end
    if isa(element, Int) && >(element, lastindex(Mendeleev.elements))
        throw(ArgumentError("Provided atomic number is > $(lastindex(Mendeleev.elements))"))
    end
    radius = getindex(Mendeleev.elements, element).ionic_radii(;
        coordination = coordination,
        charge = charge
    )[1].ionic_radius
    if value_only
        return radius.val
    else
        return radius
    end
end

"""
    atomic_mass(element::Union{String,Symbol,Int}; [value_only::Bool = false])

Helper function to return the unified atomic mass (u) of an element.

Wraps around Mendeleev.elements[:Symbol], see [Mendeleev.jl](https://eben60.github.io/Mendeleev.jl/)

# Examples

```julia-repl
julia> atomic_mass(:U)
238.02891 u

julia> atomic_mass(\"U\")
238.02891 u
```
"""
function atomic_mass(element::Union{String,Symbol,Int}; value_only::Bool = false)
    if isa(element, AbstractString)
        element = Symbol(element)
    end
    if isa(element, Int) && >(element, lastindex(Mendeleev.elements))
        throw(ArgumentError("Provided atomic number is > $(lastindex(Mendeleev.elements))"))
    end
    mass = getindex(Mendeleev.elements, element).atomic_mass
    if value_only
        return mass.val
    else
        return mass
    end
end
