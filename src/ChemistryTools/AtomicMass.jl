# SPDX-FileCopyrightText: 2024 Jarred Lloyd (https://github.com/jarredclloyd)
# SPDX-License-Identifier: MIT
# This file contains a dictionary of the IUPAC aggregate atomic masses for each element and a function to help the user retrieve them.

export get_atomicmass

"""
    get_atomicmass(element)

    Returns the atomic mass of an element from the PeridocTables.jl package.

# Argument Notes
    -  The only required argument `element` can be entered as a (case insensitive) string, character, symbol or integer representing the atomic number, symbol or name of an element.

```julia-repl
julia> get_atomicmass('H')
1.008 u
julia> get_atomicmass("He")
4.002602 u
julia> get_atomicmass("Lithium")
6.94 u
julia> get_atomicmass("beryllium")
9.0121831 u
julia> get_atomicmass(:B)
10.81 u
julia> get_atomicmass(6)
12.001 u
```
"""
function get_atomicmass(element::Union{T} where T <: Union{AbstractString, AbstractChar, Symbol, Integer})
    # Convert the input to uppercase to handle both symbol and name inputs uniformly.
    if typeof(element) <: Symbol
    elseif (typeof(element) <: Integer)
    elseif length(element) ≤ 2
        element = Symbol(element)
    end
    if typeof(element) <: Integer && element in keys(PeriodicTable.elements.bynumber)
        return elements[element].atomic_mass
    elseif typeof(element) == Symbol && element in keys(PeriodicTable.elements.bysymbol)
        return elements[element].atomic_mass
    elseif typeof(element) <: AbstractString && lowercase(element) in keys(PeriodicTable.elements.byname)
        return elements[lowercase(element)].atomic_mass
    else
        throw(ArgumentError("Element $element could not be found by symbol, name, or atomic number in the `PeriodicTable.elements` dictionary."))
    end
end
