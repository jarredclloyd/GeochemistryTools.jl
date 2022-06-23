#= 
MOLECULARWT - Computes the molecular weight for a given chemical formula.

weight = MOLECULARWT(F) where F is a string giving a chemical formula unit (e.g. Fe2O3). Correct capitalization is essential to correctly interpreting the formula. The molecular weight, weight, is returned in g/mol.

If F is an array, each element of F is a separate formula string.
=#

# Dictionaries, mass in g/mol
element_symbol_to_mass = Dict([('H', 1.00794), ("He", 4.002602), ("Li", 6.941), ("Be", 9.012182), ('B', 10.811), ('C', 12.0107), ('N', 14.0067), ('O', 15.9994), ('F', 18.9984032), ("Ne", 20.1797), ("Na", 22.98976928), ("Mg", 24.305), ("Al", 26.9815386), ("Si", 28.0855), ('P', 30.973762), ('S', 32.065), ("Cl", 35.453), ("Ar", 39.948), ('K', 39.0983), ("Ca", 40.078), ("Sc", 44.955912), ("Ti", 47.867), ('V', 50.9415), ("Cr", 51.9961), ("Mn", 54.938045), ("Fe", 55.845), ("Co", 58.933195), ("Ni", 58.6934), ("Cu", 63.546), ("Zn", 65.409), ("Ga", 69.723), ("Ge", 72.64), ("As", 74.9216), ("Se", 78.96), ("Br", 79.904), ("Kr", 83.798), ("Rb", 85.4678), ("Sr", 87.62), ('Y', 88.90585), ("Zr", 91.224), ("Nb", 92.90638), ("Mo", 95.94), ("Ru", 101.07), ("Rh", 102.9055), ("Pd", 106.42), ("Ag", 107.8682), ("Cd", 112.411), ("In", 114.818), ("Sn", 118.71), ("Sb", 121.76), ("Te", 127.6), ('I', 126.90447), ("Xe", 131.293), ("Cs", 132.9054519), ("Ba", 137.327), ("La", 138.90547), ("Ce", 140.116), ("Pr", 140.90765), ("Nd", 144.242), ("Sm", 150.36), ("Eu", 151.964), ("Gd", 157.25), ("Tb", 158.92535), ("Dy", 162.5), ("Ho", 164.93032), ("Er", 167.259), ("Tm", 168.93421), ("Yb", 173.04), ("Lu", 174.967), ("Hf", 178.49), ("Ta", 180.94788), ("W", 183.84), ("Re", 186.207), ("Os", 190.23), ("Ir", 192.217), ("Pt", 195.084), ("Au", 196.966569), ("Hg", 200.59), ("Tl", 204.3833), ("Pb", 207.2), ("Bi", 208.9804), ("Th", 232.03806), ("Pa", 231.03588), ('U', 238.0289)])

#Primary function
function molecularwt(formula)
    if isa(formula, Array)
        weight = zeros(size(formula))
        formulas_count = length(formula)
        for i in 1:formulas_count
            weight[i] = computemass(formula[i])
        end
    else
        weight = computemass(formula)
    end
    if isa(formula, Array)
        for i in 1:formulas_count
            println("Molecular weight of input formula {" * string(formula[i]) * "} is " * string(weight[i]) * " g/mol")
        end
    else
    println("Molecular weight of input formula {" * string(formula) * "} is " * string(weight) * " g/mol")
    end
end

#Compute formula mass
function computemass(formula)
    weight = 0
    position_index = 1
    formula_length = length(formula)
    while position_index <= formula_length
        mass, position_index = computeionmass(formula, position_index)
        weight = weight + mass
    end
    return weight
end

#Find ionic masses
function computeionmass(formula, left_index)
    formula_length = length(formula)
    right_index = left_index + 1
    if right_index < formula_length && formula[right_index] === '(' #to catch ")(' in formulas
        left_index = left_index + 1
        right_index = left_index + 1 
    end
    if formula[left_index] === '('  #to catch parenthesised ions
        left_index = left_index + 1
        if findnext(')', formula, left_index) !== nothing
            right_index = findnext(')', formula, left_index)
            parentheses_end = right_index
            if right_index - left_index === 1
                mass = element_symbol_to_mass[formula[left_index]]
                println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(mass) * " g/mol")
                left_index = left_index + 2
            elseif right_index - left_index === 2    #e.g., (OH) (Mg) (B2)
                right_index = right_index - 1
                if formula[right_index] == uppercase(formula[right_index]) && tryparse(Int, string(formula[right_index])) === nothing
                    ion_mass = element_symbol_to_mass[formula[left_index]]
                    println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(ion_mass) * " g/mol")
                    ion_mass = ion_mass + element_symbol_to_mass[formula[right_index]]
                    println("Last calculation = " * string(formula[right_index]) * " @ total weight of " * string(ion_mass) * " g/mol")
                    left_index = left_index + 3
                elseif formula[right_index] === uppercase(formula[right_index]) && tryparse(Int, string(formula[right_index])) !== nothing
                    ion_mass = element_symbol_to_mass[formula[left_index]]
                    n_moles = parse(Int, formula[right_index])
                    ion_mass = n_moles * ion_mass
                    println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(ion_mass) * " g/mol")
                    left_index = left_index + 3
                else
                    ion_mass = element_symbol_to_mass[formula[left_index:right_index]]
                    println("Last calculation = " * string(formula[left_index:right_index]) * " @ total weight of " * string(ion_mass) * " g/mol")
                    left_index = left_index + 3
                end
                n_moles, position_index = findnmoles(formula, left_index - 1)
                mass = ion_mass * n_moles
            elseif right_index - left_index > 2     #e.g., (Mg2Si3O10OH)
                right_index = left_index + 1
                mass = 0   #reset mass before cumulative calculation through parenthesised ions
                while right_index <= parentheses_end
                    if formula[right_index] === ')' && tryparse(Int, string(formula[left_index])) === nothing   #e.g., H)
                        ion_mass = element_symbol_to_mass[formula[left_index]]
                        println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(ion_mass) * " g/mol")
                        left_index = right_index + 1
                        right_index = left_index + 1
                        position_index = left_index
                    elseif tryparse(Int, string(formula[right_index])) !== nothing && formula[left_index] === uppercase(formula[left_index])   #e.g., O10, O3
                        ion_mass = element_symbol_to_mass[formula[left_index]]
                        n_moles, position_index = findnmoles(formula, right_index)
                        ion_mass = n_moles * ion_mass
                        println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(ion_mass) * " g/mol")
                        left_index = position_index
                        right_index = left_index + 1
                    elseif tryparse(Int, string(formula[right_index])) === nothing && formula[right_index] === lowercase(formula[right_index])  #e.g., Mg
                        ion_mass = element_symbol_to_mass[formula[left_index:right_index]]
                        n_moles, position_index = findnmoles(formula, right_index + 1)
                        ion_mass = n_moles * ion_mass
                        println("Last calculation = " * string(formula[left_index:right_index]) * " @ total weight of " * string(ion_mass) * " g/mol")
                        left_index = position_index
                        right_index = left_index + 1
                    elseif tryparse(Int, string(formula[right_index])) === nothing && formula[right_index] === uppercase(formula[right_index]) #e.g., O in OH
                        ion_mass = element_symbol_to_mass[formula[left_index]]
                        println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(ion_mass) * " g/mol")
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
                error("ERROR: Formula is malformed. Length between parentheses at positions" * string(left_index) * " and " * string(right_index) * " is 1 or 0.")
            end
        else
            error("ERROR: Missing closing parentheses after string position" * string(left_index))
        end
        n_moles, position_index = findnmoles(formula, left_index)
        mass = mass * n_moles #final mass within parentheses multiplied by n moles of parenthesised ions e.g., )2
        left_index = position_index
        println("Total mass of preceding parenthesised ions in formula is " *string(mass) *"g/mol")
    elseif right_index <= formula_length
        if tryparse(Int, string(formula[right_index])) !== nothing && formula[left_index] === uppercase(formula[left_index])    #e.g., B2
            ion_mass = element_symbol_to_mass[formula[left_index]]
            n_moles, position_index = findnmoles(formula, right_index)
            mass = ion_mass * n_moles
            println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(mass) * " g/mol")
        elseif tryparse(Int, string(formula[right_index])) === nothing && formula[right_index] === lowercase(formula[right_index])     #e.g., Mg
            ion_mass = element_symbol_to_mass[formula[left_index:right_index]]
            n_moles, position_index = findnmoles(formula, right_index + 1)
            mass = ion_mass * n_moles
            println("Last calculation = " * string(formula[left_index:right_index]) * " @ total weight of " * string(mass) * " g/mol")
        elseif formula[right_index] == uppercase(formula[right_index])      #e.g., O in MgOH
            ion_mass = element_symbol_to_mass[formula[left_index]]
            mass = ion_mass
            position_index = right_index
            println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(mass) * " g/mol")
        end
    elseif right_index > formula_length && formula[left_index] === uppercase(formula[left_index])   #e.g., O  in MgO
        ion_mass = element_symbol_to_mass[formula[left_index]]
        mass = ion_mass
        position_index = right_index
        println("Last calculation = " * string(formula[left_index]) * " @ total weight of " * string(mass) * " g/mol")
    end
    return mass, position_index
end

#Find number of moles for ion
function findnmoles(formula, mole_index)
    formula_length = length(formula)
    if mole_index + 2 <= formula_length && formula[mole_index+2] === '.' #find decimals >10.00
        if mole_index + 4 <= formula_length && tryparse(Int, string(formula[mole_index+4])) !== nothing
            n_moles = parse(Float64, formula[mole_index:mole_index+4])
            position_index = mole_index + 5
        elseif mole_index + 3 <= formula_length && tryparse(Int, string(formula[mole_index+3])) !== nothing
            n_moles = parse(Float64, formula[mole_index:mole_index+2])
            position_index = mole_index + 4
        end
    elseif mole_index + 1 <= formula_length && formula[mole_index+1] !== '.' #find Int
        if tryparse(Int, string(formula[mole_index+1])) !== nothing && tryparse(Int, string(formula[mole_index])) !== nothing  #e.g., Si10
            n_moles = parse(Int, formula[mole_index:mole_index+1])
            position_index = mole_index + 2
        elseif tryparse(Int, string(formula[mole_index])) !== nothing  #e.g., B2
            n_moles = parse(Int, formula[mole_index])
            position_index = mole_index + 1
        else
            n_moles = 1
            position_index = mole_index
        end
    elseif mole_index + 1 <= formula_length && formula[mole_index+1] === '.' #find decimals <10.00
        if mole_index + 3 <= formula_length && tryparse(Int, string(formula[mole_index+3])) !== nothing
            n_moles = parse(Float64, formula[mole_index:mole_index+3])
            position_index = mole_index + 4
        elseif mole_index + 2 <= formula_length && tryparse(Int, string(formula[mole_index+2])) !== nothing
            n_moles = parse(Float64, formula[mole_index:mole_index+2])
            position_index = mole_index + 3
        end
    elseif mole_index <= formula_length && tryparse(Int, string(formula[mole_index])) !== nothing
        n_moles = parse(Int, string(formula[mole_index]))
        position_index = mole_index + 1
    else
        n_moles = 1
        position_index = mole_index
    end
    if position_index <= formula_length && formula[position_index] === ')'
        position_index = position_index + 1
    end
    return n_moles, position_index
end