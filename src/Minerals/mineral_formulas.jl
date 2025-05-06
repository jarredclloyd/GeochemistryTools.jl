#= Preamble
This file contains functions for calculating mineral formulas from ppm or wt% data
=#
export formula_mica

"""
    formula_mica(data, units; [normalising_O::Int=11])

Calculates mica formula based on input data and declared units.

# Description
Data is a data frame containing sample name and elemental data (either as element or oxide) in a consistent unit
(either ppm OR wt%O).

Units is the units the data is in, a string of value "ppm" OR "wt%O"

# Example
```julia-repl
julia> formula_mica(df, "ppm")
```
"""
function formula_mica(data::AbstractDataFrame, units::AbstractString; normalising_O::Int=11)
    I_site = [:NH4, :Na, :K, :Ca, :Rb, :Cs, :Ba]
    M_site = [:Li, :Mg, :Ti, :V, :Cr, :MnO, :Mn2O3, :Zn, :FeO, :Fe2O3, :Al]
    T_site = [:Fe2O3, :Al, :Be, :B, :Si]
    A_site = [:F, :Cl, :OH, :S]
    column_names = [:NH4, :Na, :K, :Ca, :Rb, :Cs, :Ba, :Li, :Mg, :Ti, :V, :Cr, :MnO, :Mn2O3, :Zn, :FeO, :Fe2O3, :Al, :Be, :B, :Si, :F, :Cl, :OH, :S]
    workingdata = _find_columns_mica(data)
    workingvector = collect(values(workingdata[1, Not(:Sample)]))

    if lowercase(units) == "ppm"
        oxide_conversion_factors = [
            cation_to_oxide(1; cation="NH4", cation_multiplicity = 2, oxide = "(NH4)2O"),
            element_to_oxide(1; element=:Na, element_multiplicity=2, oxide = "Na2O").val,
            element_to_oxide(1; element=:K, element_multiplicity=2, oxide = "K2O").val,
            element_to_oxide(1; element=:Ca, element_multiplicity=1, oxide = "CaO").val,
            element_to_oxide(1; element=:Rb, element_multiplicity=2, oxide = "Rb2O").val,
            element_to_oxide(1; element=:Cs, element_multiplicity=2, oxide = "Cs2O").val,
            element_to_oxide(1; element=:Ba, element_multiplicity= 1, oxide = "BaO").val,
            element_to_oxide(1; element=:Li, element_multiplicity= 2, oxide = "Li2O").val,
            element_to_oxide(1; element=:Mg, element_multiplicity= 1, oxide = "MgO").val,
            element_to_oxide(1; element=:Ti, element_multiplicity= 1, oxide = "TiO2").val,
            element_to_oxide(1; element=:V, element_multiplicity= 1, oxide = "VO2").val,
            element_to_oxide(1; element=:Cr, element_multiplicity= 2, oxide = "Cr2O3").val,
            element_to_oxide(1; element=:Mn, element_multiplicity= 1, oxide = "MnO").val,
            element_to_oxide(1; element=:Mn, element_multiplicity= 2, oxide ="Mn2O3").val,
            element_to_oxide(1; element=:Zn, element_multiplicity= 1, oxide = "ZnO").val,
            element_to_oxide(1; element=:Fe, element_multiplicity=1, oxide = "FeO").val,
            element_to_oxide(1; element=:Fe, element_multiplicity=2, oxide = "Fe2O3").val,
            element_to_oxide(1; element=:Al, element_multiplicity= 2, oxide = "Al2O3").val,
            element_to_oxide(1; element=:Be, element_multiplicity= 1, oxide = "BeO").val,
            element_to_oxide(1; element=:B, element_multiplicity= 2, oxide = "B2O3").val,
            element_to_oxide(1; element=:Si, element_multiplicity= 1, oxide = "SiO2").val,
            element_to_oxide(1; element=:F, element_multiplicity= 1, oxide = "F").val,
            element_to_oxide(1; element=:Cl, element_multiplicity= 1, oxide = "Cl").val,
            1, #OH
            element_to_oxide(1; element=:S, element_multiplicity= 1, oxide = "S").val
        ]
        weight_percent = (workingvector ./ 1e4) .* oxide_conversion_factors
    elseif lowercase(units) == "wt%o"
        weight_percent = workingvector
    end

    molecular_weights = [
        molecular_mass("(NH4)2O"; verbose = false),
        molecular_mass("Na2O"; verbose = false),
        molecular_mass("K2O"; verbose = false),
        molecular_mass("CaO"; verbose = false),
        molecular_mass("Rb2O"; verbose = false),
        molecular_mass("Cs2O"; verbose = false),
        molecular_mass("BaO"; verbose = false),
        molecular_mass("Li2O"; verbose = false),
        molecular_mass("MgO"; verbose = false),
        molecular_mass("TiO2"; verbose = false),
        molecular_mass("VO2"; verbose = false),
        molecular_mass("Cr2O3"; verbose = false),
        molecular_mass("MnO"; verbose = false),
        molecular_mass("Mn2O3"; verbose = false),
        molecular_mass("ZnO"; verbose = false),
        molecular_mass("FeO"; verbose = false),
        molecular_mass("Fe2O3"; verbose = false),
        molecular_mass("Al2O3"; verbose = false),
        molecular_mass("BeO"; verbose = false),
        molecular_mass("B2O3"; verbose = false),
        molecular_mass("SiO2"; verbose = false),
        get_atomicmass(:F).val,
        get_atomicmass(:Cl).val,
        molecular_mass("OH"; verbose = false),
        get_atomicmass(:S).val
    ]

    moles_compound = weight_percent ./ molecular_weights
    oxy_per_elem = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 1, 3, 1, 1, 3, 3, 1, 3, 2, 0.5, 0.5, 0.5, 0] # need to work out S properly, this will currently zero it.
    moles_oxygen = moles_compound .* oxy_per_elem
    oxygen_sum = sum(moles_oxygen)
    normalising_factor = normalising_O / oxygen_sum
    normalised_oxygen_moles = moles_oxygen .* normalising_factor

    # iterate OH
    ϵ = 1
    OH_initial = normalised_oxygen_moles[24]

    if OH_initial == 0
        iteration = 1
        while ϵ > 1e-4 && iteration <= 1000
            OH = 2 - normalised_oxygen_moles[22] - normalised_oxygen_moles[23] - normalised_oxygen_moles[25]
            println("H iteration: $OH")
            moles_oxygen[24] = OH
            oxygen_sum = sum(moles_oxygen)
            normalising_factor = normalising_O / oxygen_sum
            normalised_oxygen_moles = moles_oxygen .* normalising_factor
            iteration += 1
            ϵ = abs(OH_initial - OH)
            println("epsilon: $ϵ")
            OH_initial = OH
        end
    end
    element_oxy_ratio = [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0.667, 2, 0.667, 2, 2, 0.667, 0.667, 2, 0.667, 1, 0.5, 0.5, 1, 1]
    atoms_per_formula_unit = normalised_oxygen_moles .* element_oxy_ratio
    println("Cation sum: $(sum(atoms_per_formula_unit[1:21]))")
    println("Total ions: $(sum(atoms_per_formula_unit))")
    return [column_names, round.(atoms_per_formula_unit,digits=3)]
end

function _find_columns_mica(data)
    workingdata = DataFrame()
    if in("sample", lowercase.(names(data)))
        workingdata.Sample = data[:, findfirst(lowercase.(names(data)) .== "sample")]
    else
        workingdata.Sample .= "sample"
    end
    if in("NH4",names(data))
        workingdata.NH4 = data[:, :NH4]
    else
        workingdata.NH4 .= 0
    end
    if in("Na",names(data))
        workingdata.Na = data[:, :Na]
    elseif in("Na2O",names(data))
        workingdata.Na = data[:, :Na2O]
    else
        workingdata.Na .= 0
    end
    if in("K",names(data))
        workingdata.K = data[:, :K]
    elseif in("K2O",names(data))
        workingdata.K = data[:, :K2O]
    else
        workingdata.K .= 0
    end
    if in("Ca",names(data))
        workingdata.Ca = data[:, :Ca]
    elseif in("CaO",names(data))
        workingdata.Ca = data[:, :CaO]
    else
        workingdata.Ca .= 0
    end
    if in("Rb",names(data))
        workingdata.Rb = data[:, :Rb]
    elseif in("Rb2O",names(data))
        workingdata.Rb = data[:, :Rb2O]
    else
        workingdata.Rb .= 0
    end
    if in("Cs",names(data))
        workingdata.Cs = data[:, :Cs]
    elseif in("Cs2O",names(data))
        workingdata.Cs = data[:, :Cs2O]
    else
        workingdata.Cs .= 0
    end
    if in("Ba",names(data))
        workingdata.Ba = data[:, :Ba]
    elseif in("BaO",names(data))
        workingdata.Ba = data[:, :BaO]
    else
        workingdata.Ba .= 0
    end
    if in("Li",names(data))
        workingdata.Li = data[:, :Li]
    elseif in("Li2O",names(data))
        workingdata.Li = data[:, :Li2O]
    else
        workingdata.Li .= 0
    end
    if in("Mg",names(data))
        workingdata.Mg = data[:, :Mg]
    elseif in("MgO",names(data))
        workingdata.Mg = data[:, :MgO]
    else
        workingdata.Mg .= 0
    end
    if in("Ti",names(data))
        workingdata.Ti = data[:, :Ti]
    elseif in("TiO2",names(data))
        workingdata.Ti = data[:, :TiO2]
    else
        workingdata.Ti .= 0
    end
    if in("V",names(data))
        workingdata.V = data[:, :V]
    elseif in("VO2",names(data))
        workingdata.V = data[:, :VO2]
    else
        workingdata.V .= 0
    end
    if in("Cr",names(data))
        workingdata.Cr = data[:, :Cr]
    elseif in("Cr2O3",names(data))
        workingdata.Cr = data[:, :Cr2O3]
    else
        workingdata.Cr .= 0
    end
    if in("MnO",names(data))
        workingdata.MnO = data[:, :MnO]
    else
        workingdata.MnO .= 0
    end
    if in("Mn2O3",names(data))
        workingdata.Mn2O3 = data[:, :Mn2O3]
    else
        workingdata.Mn2O3 .= 0
    end
    if in("Mn",names(data))
        workingdata.MnO .= 0
    end
    if in("Zn",names(data))
        workingdata.Zn = data[:, :Zn]
    elseif in("ZnO",names(data))
        workingdata.Zn = data[:, :ZnO]
    else
        workingdata.Zn .= 0
    end
    if in("FeO",names(data))
        workingdata.FeO = data[:, :FeO]
    else
        workingdata.FeO .= 0
    end
    if in("Fe2O3",names(data))
        workingdata.Fe2O3 = data[:, :Fe2O3]
    else
        workingdata.Fe2O3 .= 0
    end
    if in("Fe",names(data))
        workingdata.FeO = data[:, :Fe]
    end
    if in("Al",names(data))
        workingdata.Al = data[:, :Al]
    elseif in("Al2O3",names(data))
        workingdata.Al = data[:, :Al2O3]
    else
        workingdata.Al .= 0
    end
    if in("Be",names(data))
        workingdata.Be = data[:, :Be]
    elseif in("BeO",names(data))
        workingdata.Be = data[:, :BeO]
    else
        workingdata.Be .= 0
    end
    if in("B",names(data))
        workingdata.B = data[:, :B]
    elseif in("B2O3",names(data))
        workingdata.B = data[:, :B2O3]
    else
        workingdata.B .= 0
    end
    if in("Si",names(data))
        workingdata.Si = data[:, :Si]
    elseif in("SiO2",names(data))
        workingdata.Si = data[:, :SiO2]
    else
        workingdata.Si .= 0
    end
    if in("F",names(data))
        workingdata.F = data[:, :F]
    else
        workingdata.F .= 0
    end
    if in("Cl",names(data))
        workingdata.Cl = data[:, :Cl]
    else
        workingdata.Cl .= 0
    end
    if in("OH",names(data))
        workingdata.OH = data[:, :OH]
    else
        workingdata.OH .= 0
    end
    if in("S",names(data))
        workingdata.S = data[:, :S]
    else
        workingdata.S .= 0
    end
    return workingdata
end

SSP18 = DataFrame([:SiO2 =>34.79,	:TiO2=>3.26,:Al2O3=>18.82,	:FeO=>21.39,	:MnO=>0.51,	:MgO=>7.62,	:CaO=>0,	:Na2O=>0.12,	:K2O=>9.66,	:BaO=>0.14,	:F=>0.17,	:Cl=>0.05])
