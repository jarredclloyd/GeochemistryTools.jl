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
    println(typeof(formula))
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
function computemoleculemass(formula)
    moleculemass = 0
    position_index = firstindex(formula)
    formula_length = length(formula)
    if length(findall('(', formula)) != length(findall(')', formula))
        error("Error: Mismatched parentheses `(` & `)` in formula")
    end
    if length(findall('[', formula)) != length(findall(']', formula))
        error("Error: Mismatched square brackets `[` & `]` in formula")
    end
    while position_index ≤ formula_length
        ionmass, position_index = computeionmass(formula, position_index, formula_length)
        moleculemass = moleculemass + ionmass
    end
    return moleculemass
end

#Find ions and compute their masses
function computeionmass(formula, left_index, formula_length)
    if left_index < lastindex(formula)
        right_index = nextind(formula, left_index)
    else
        right_index = lastindex(formula)
    end
    if formula[left_index] == '['
        sqrbrac_open = formula[left_index]
        sqrbrac_close = formula[findnext('[', formula, left_index)]
        left_index = nextind(formula, sqrbrac_open)
        right_index = nextind(formula, sqrbrac_open)
    end
    if formula[left_index] == '('
        paren_open = formula[left_index]
        paren_close = formula[findnext(')', formula, left_index)]
        left_index = nextind(formula, paren_open)
        right_index = nextind(formula, paren_open)
    end
    if right_index <= formula_length &&
    occursin(r"[A-Z][a-z]", formula[left_index:right_index]) &&
    in(Symbol(formula[left_index:right_index]), keys(PeriodicTable.elements.bysymbol))
        ionicmass = get_atomicmass(formula[left_index:right_index]).val
    elseif occursin(r"[A-Z]", string(formula[left_index])) &&
    in(Symbol(formula[left_index]), keys(PeriodicTable.elements.bysymbol))
       ionicmass = get_atomicmass(formula[left_index]).val
    end
    if checkbounds(Bool, formula, nextind(formula, right_index))
        left_index = nextind(formula, right_index)
        if left_index < lastindex(formula)
            if isnumeric(formula[left_index])
                right_index = findnext(r"[A-Z()\[\]]",formula, left_index)[1]
                n_mole_end_ind = right_index - 1
                if n_mole_end_ind == left_index
                    ionicmass *= tryparse(Int64, formula[left_index])
                else
                    ionicmass *= tryparse(Float64, formula[left_index:n_mole_end_ind])
                end
                left_index = right_index
                right_index = nextind(formula, left_index)
            end
            if left_index < paren_close
                cummulativemass = ionicmass
                while left_index < paren_close
                    ionicmass = computeionmass(formula[left_index:right_index], left_index, length(formula[left_index:paren_close]))
                    cummulativemass += ionicmass
                    left_index = right_index
                    right_index = nextind(formula, left_index)
                end
            end
        end
    end
    position_index = right_index+1
    return ionicmass, position_index
end

#Find number of moles for ion
function findnmoles(formula, mole_index)
    formula_length = length(formula)
    #if findnext('.', formula, mole_index) != nothing
    #    decimal_index = findnext('.', formula, mole_index)
    if mole_index + 2 ≤ formula_length && formula[mole_index+2] == '.' #find decimals >10.00
        if mole_index + 4 ≤ formula_length && isdigit(formula[mole_index+4]) == true
            n_moles = parse(Float64, formula[mole_index:mole_index+4])
            position_index = mole_index + 5
        elseif mole_index + 3 ≤ formula_length && isdigit(formula[mole_index+3]) == true
            n_moles = parse(Float64, formula[mole_index:mole_index+2])
            position_index = mole_index + 4
        end
    elseif mole_index + 1 ≤ formula_length && formula[mole_index+1] != '.' #find Int
        if isdigit(formula[mole_index+1]) == true && isdigit(formula[mole_index]) == true  #e.g., Si10
            n_moles = parse(Int, formula[mole_index:mole_index+1])
            position_index = mole_index + 2
        elseif isdigit(formula[mole_index]) == true  #e.g., B2
            n_moles = parse(Int, formula[mole_index])
            position_index = mole_index + 1
        else
            n_moles = 1
            position_index = mole_index
        end
    elseif mole_index + 1 ≤ formula_length && formula[mole_index+1] == '.' #find decimals <10.00
        if mole_index + 3 ≤ formula_length && isdigit(formula[mole_index+3]) == true
            n_moles = parse(Float64, formula[mole_index:mole_index+3])
            position_index = mole_index + 4
        elseif mole_index + 2 ≤ formula_length && isdigit(formula[mole_index+2]) == true
            n_moles = parse(Float64, formula[mole_index:mole_index+2])
            position_index = mole_index + 3
        end
    elseif mole_index ≤ formula_length && isdigit(formula[mole_index]) == true
        n_moles = parse(Int, string(formula[mole_index]))
        position_index = mole_index + 1
    else
        n_moles = 1
        position_index = mole_index
    end
    if position_index ≤ formula_length && formula[position_index] == ')'
        position_index = position_index + 1
    end
    return n_moles, position_index
end
