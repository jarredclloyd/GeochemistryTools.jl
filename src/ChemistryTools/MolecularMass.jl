#= Preamble
This file contains various decay system equations for calculation of ratios and ages (U-Pb).
=#

export molecular_mass

"""
    molecular_mass(formula; verbose::Bool=true)

Calculates molar weight of input formula.

Multiple formulas can be passed as vector of strings.
Use only standard characters and do not superscript mole numbers. Additional ions after a • should be added in
parentheses if there are multiple moles of it, e.g., •2(H20) should be entered as (H20)2.

# Keywords
    `-verbose::Bool=true` set to false to suppress detailed output to REPL.

```julia-repl
julia> molecular_mass("Na0.84Ca0.64K0.44Ba0.02Sr0.01Ca0.64Mn0.05Zn0.01Nb2.7Ti1.28Fe0.07Si7.97Al0.06O24.09OH1.23O2.75(H2O)6.92")
Molecular weight of input formula {Na0.84Ca0.64K0.44Ba0.02Sr0.01Ca0.64Mn0.05Zn0.01Nb2.7Ti1.28Fe0.07Si7.97Al0.06O24.09OH1.23O2.75(H2O)6.92} is:
1207.6534935612 gmol⁻¹.

julia> molecular_mass(["Al2O3", "K2O"])
Molecular weight of input formula {K2O} is:
94.196 gmol⁻¹.
Molecular weight of input formula {Al2O3} is:
101.9612772 gmol⁻¹.
```
"""
function molecular_mass(formula::Union{String, Vector{String}};
    verbose::Bool=true)
    if isa(formula, Vector{String})
        mass = zeros(size(formula))
        nformulas = length(formula)
        for i in 1:nformulas
            mass[i] = computemoleculemass(formula[i])
            if verbose == true
                printstyled("Molecular mass of input formula {" * string(formula[i]) * "} is:\n" * string(mass[i]) *
                        " gmol⁻¹.\n", bold=true)
            end
        end
        return mass
    else
        mass = computemoleculemass(formula)
        if verbose == true
        printstyled("Molecular mass of input formula {" * string(formula) * "} is:\n" * string(mass) *
                    " gmol⁻¹.", bold=true)
        end
        return mass
    end
end

#Compute total formula mass
function computemoleculemass(formula::AbstractString)
    moleculemass = 0
    position_index = firstindex(formula)
    formula_lastind = lastindex(formula)
    if length(findall('(', formula)) != length(findall(')', formula))
        error("Error: Mismatched parentheses `(` & `)` in formula")
    end
    if length(findall('[', formula)) != length(findall(']', formula))
        error("Error: Mismatched square brackets `[` & `]` in formula")
    end
    while position_index ≤ formula_lastind
        ionmass, position_index = computeionmass(formula, position_index, formula_lastind)
        moleculemass = moleculemass + ionmass
    end
    return moleculemass
end

#Find ions and compute their masses
function computeionmass(formula::AbstractString, left_index::Integer, formula_lastind::Integer)
    parenmass = nothing
    sqrbracmass = nothing
    if left_index < lastindex(formula)
        right_index = nextind(formula, left_index)
    else
        right_index = lastindex(formula)
    end
    # find parentheses and brackets
    if formula[left_index] == '['
        sqrbrac_close = findnext('[', formula, left_index)[1]
        sqrbracmass, left_index = _sqrbracsum(formula, left_index+1, sqrbrac_close,formula_lastind)
        left_index = nextind(formula, left_index)
    end
    if formula[left_index] == '('
        paren_close = findnext(')', formula, left_index)[1]
        parenmass, left_index =_parensum(formula, left_index+1, paren_close,formula_lastind)
        left_index = nextind(formula, left_index)
    end
    if isnothing(sqrbracmass) && isnothing(parenmass)
        element = _get_element!(formula, left_index, formula_lastind)
        ionicmass = get_atomicmass(element).val
        left_index = left_index + length(element)
    elseif !isnothing(parenmass)
        ionicmass = parenmass
    elseif !isnothing(sqrbracmass)
        ionicmass = sqrbracmass
    end
    if left_index ≤ formula_lastind
        if isnumeric(formula[left_index])
            nmoles, left_index = _get_nmoles(formula, left_index, formula_lastind)
            ionicmass *= nmoles
        end
    end
    position_index = left_index
    return ionicmass, position_index
end


function _get_element!(formula::AbstractString, elem_start_ind::Integer, formula_lastind::Integer)
    element = nothing
    if occursin(r"[()\[\]]", string(formula[elem_start_ind]))
    else
        if elem_start_ind < formula_lastind
            elem_end_ind = nextind(formula, elem_start_ind)
        else
            elem_end_ind = elem_start_ind
        end
        if occursin(r"[0-9.()\[\]]", string(formula[elem_end_ind]))
            elem_end_ind = elem_start_ind
        end
        if elem_start_ind < formula_lastind && occursin(r"[A-Z][a-z]", formula[elem_start_ind:elem_end_ind]) &&
            in(Symbol(formula[elem_start_ind:elem_end_ind]), keys(PeriodicTable.elements.bysymbol))
            return element = formula[elem_start_ind:elem_end_ind]
        elseif elem_start_ind ≤ formula_lastind &&
            occursin(r"[A-Z]", string(formula[elem_start_ind])) &&
            in(Symbol(formula[elem_start_ind]), keys(PeriodicTable.elements.bysymbol))
            return element = formula[elem_start_ind]
        end
        if isnothing(element)
                error("Invalid element symbol ($(formula[elem_start_ind:elem_end_ind])) found in formula at position $elem_start_ind")
        else
        return element
        end
    end
end

function _get_nmoles(formula::AbstractString, nmole_start_ind::Integer, formula_lastind::Integer)
    if nmole_start_ind == formula_lastind
        nmoles = tryparse(Int, string(formula[nmole_start_ind]))
        nmole_end_ind = nextind(formula, nmole_start_ind)+1
    else
        nmole_end_ind = findnext(r"[A-Za-z()\[\]]",formula, nmole_start_ind)
        if isnothing(nmole_end_ind)
            nmole_end_ind = formula_lastind
            nmoles = tryparse(Float64, string(formula[nmole_start_ind:nmole_end_ind]))
            nmole_end_ind += 1
        else
            nmole_end_ind = nmole_end_ind[1]
            nmoles = tryparse(Float64, string(formula[nmole_start_ind:nmole_end_ind-1]))
        end
    end
    return nmoles, nmole_end_ind
end


function _parensum(formula::AbstractString, left_index::Integer, paren_close::Integer, formula_lastind::Integer)
    cumulativemass = 0.0
    while left_index < paren_close
        ionmass, left_index = computeionmass(formula, left_index, formula_lastind::Integer)
        cumulativemass += ionmass
    end
    return cumulativemass, left_index
end

function _sqrbracsum(formula::AbstractString, left_index::Integer, sqrbrac_close::Integer, formula_lastind::Integer)
    cumulativemass = 0.0
    while left_index < sqrbrac_close
        ionmass, left_index = computeionmass(formula, left_index, formula_lastind::Integer)
        cumulativemass += ionmass
    end
    return cumulativemass, left_index
end
