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
        weight = zeros(size(formula))
        nformulas = length(formula)
        for i in 1:nformulas
            weight[i] = computemass(formula[i])
            if verbose == true
                printstyled("Molecular weight of input formula {" * string(formula[i]) * "} is:\n" * string(weight[i]) *
                        " gmol⁻¹.\n", bold=true)
            end
        end
        return weight
    else
        weight = computemass(formula)
        if verbose == true
        printstyled("Molecular weight of input formula {" * string(formula) * "} is:\n" * string(weight) *
                    " gmol⁻¹.", bold=true)
        end
        return weight
    end
end

#Compute total formula mass
function computemass(formula)
    weight = 0
    position_index = 1
    formula_length = length(formula)
    while position_index ≤ formula_length
        mass, position_index = computeionmass(formula, position_index)
        weight = weight + mass
    end
    return weight
end

#Find ions and compute their masses
function computeionmass(formula, left_index)
    formula_length = length(formula)
    right_index = left_index + 1
    if right_index < formula_length && formula[right_index] == '(' #to catch ")(' in formulas
        left_index = left_index + 1
        right_index = left_index + 1
    end
    if formula[left_index] == '('  #to catch parenthesised ions
        left_index = left_index + 1
        if findnext(')', formula, left_index) != nothing
            right_index = findnext(')', formula, left_index)
            parentheses_end = right_index
            if right_index - left_index == 1
                mass = atomic_mass[formula[left_index]]
                left_index = left_index + 2
            elseif right_index - left_index == 2   #e.g., (OH) (Mg) (B2)
                right_index -= - 1
                if formula[right_index] == uppercase(formula[right_index]) && isdigit(formula[right_index]) == false
                    ion_mass = atomic_mass[formula[left_index]]
                    ion_mass = ion_mass + atomic_mass[formula[right_index]]
                    left_index = left_index + 3
                elseif formula[right_index] == uppercase(formula[right_index]) && isdigit(formula[right_index]) == true
                    ion_mass = atomic_mass[formula[left_index]]
                    n_moles = parse(Int, formula[right_index])
                    ion_mass = n_moles * ion_mass
                    left_index = left_index + 3
                else
                    ion_mass = atomic_mass[formula[left_index:right_index]]
                    left_index = left_index + 3
                end
                n_moles, position_index = findnmoles(formula, left_index - 1)
                mass = ion_mass * n_moles
            elseif right_index - left_index > 2     #e.g., (Mg2Si3O10OH)
                right_index = left_index + 1
                mass = 0   #reset mass before cumulative calculation through parenthesised ions
                while right_index ≤ parentheses_end
                    if formula[right_index] == ')' && isdigit(formula[left_index]) == false   #e.g., H)
                        ion_mass = atomic_mass[formula[left_index]]
                        left_index = right_index + 1
                        right_index = left_index + 1
                        position_index = left_index
                    elseif isdigit(formula[right_index]) == true && formula[left_index] ==
                                                                   uppercase(formula[left_index])   #e.g., O10, O3
                        ion_mass = atomic_mass[formula[left_index]]
                        n_moles, position_index = findnmoles(formula, right_index)
                        ion_mass = n_moles * ion_mass
                        left_index = position_index
                        right_index = left_index + 1
                    elseif isdigit(formula[right_index]) == false && formula[right_index] ==
                                                                    lowercase(formula[right_index])  #e.g., Mg
                        ion_mass = atomic_mass[formula[left_index:right_index]]
                        n_moles, position_index = findnmoles(formula, right_index + 1)
                        ion_mass = n_moles * ion_mass
                        left_index = position_index
                        right_index = left_index + 1
                    elseif isdigit(formula[right_index]) == false && formula[right_index] ==
                                                                    uppercase(formula[right_index]) #e.g., O in OH
                        ion_mass = atomic_mass[formula[left_index]]
                        left_index = right_index
                        right_index = left_index + 1
                        position_index = left_index
                    else
                        error("ERROR: Could not determine ion(s) in parentheses.")
                        break
                    end
                    mass = mass + ion_mass
                end
            elseif right_index - left_index < 2
                error("ERROR: Formula is malformed. Length between parentheses at positions" * string(left_index) *
                      " and " * string(right_index) * " is 1.")
            end
        else
            error("ERROR: Missing closing parentheses after string position" * string(left_index))
        end
        n_moles, position_index = findnmoles(formula, left_index)
        mass = mass * n_moles #final mass within parentheses multiplied by n moles of parenthesised ions e.g., )2
        left_index = position_index
    elseif right_index ≤ formula_length
        if isdigit(formula[right_index]) == true && formula[left_index] == uppercase(formula[left_index])    #e.g., B2
            ion_mass = atomic_mass[formula[left_index]]
            n_moles, position_index = findnmoles(formula, right_index)
            mass = ion_mass * n_moles
        elseif isdigit(formula[right_index]) == false && formula[right_index] == lowercase(formula[right_index])     #e.g., Mg
            ion_mass = atomic_mass[formula[left_index:right_index]]
            n_moles, position_index = findnmoles(formula, right_index + 1)
            mass = ion_mass * n_moles
        elseif formula[right_index] == uppercase(formula[right_index])      #e.g., O in MgOH
            ion_mass = atomic_mass[formula[left_index]]
            mass = ion_mass
            position_index = right_index
        end
    elseif right_index > formula_length && formula[left_index] == uppercase(formula[left_index])   #e.g., O  in MgO
        ion_mass = atomic_mass[formula[left_index]]
        mass = ion_mass
        position_index = right_index
    end
    return mass, position_index
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
