#= Preamble
This file contains functions for calculating mineral formulas from ppm or wt% data
=#
export formula_mica

"""
    formula_mica(data, units; [normalising_O::Int=11])

Calculates mica formula based on input data and declared units.

# Description
Data is a data frame containing sample name and elemental data (either as element or oxide) in a consistent unit
(either ppm OR wt%).

Units is the units the data is in, a string of value "ppm" OR "wt%"

# Example
```julia-repl
julia> formula_mica(df, "ppm")
```
"""
function formula_mica(data::AbstractDataFrame, units::AbstractString; normalising_O::Int=11)
    I_site = [:Na, :K, :Ca, :Rb, :Cs, :Ba]
    M_site = [:Li, :Mg, :Ti, :V, :Cr, :Mn, :Zn, :Fe]
    TM_site = [:Al]
    T_site = [:Be, :B, :Si]
    A_site = [:F, :Cl, :OH]
    column_names = [:Na, :K, :Ca, :Rb, :Cs, :Ba, :Li, :Mg, :Ti, :V, :Cr, :Mn, :Zn, :Fe, :Al, :Be, :B, :Si, :F, :Cl, :OH]
    workingdata = _find_columns_mica(data)
    workingvector = collect(values(workingdata[1, Not(:Sample)]))

    if lowercase(units) == "ppm"
        oxide_conversion = [
            1.347955633, 1.204601258, 1.399196567, 1.093596434, 1.060187345, 1.1165004, 2.152665706, 1.658259617,
            1.668477239, 1.628126104, 1.461545119, 1.291219193, 1.244707862, 1.28648939, 1.889426284, 2.775260203,
            3.220027752, 2.139327043, 1.0, 1.0, 8.936011905
        ]
        weight_percent = (workingvector ./ 1e4) .* oxide_conversion
    elseif lowercase(units) == "wt%"
        weight_percent = workingvector
    end

    molecular_weights = [
        molecular_mass("Na2O"),
        molecular_mass("K2O"),
        molecular_mass("CaO"),
        molecular_mass("Rb2O"),
        molecular_mass("Cs2O"),
        molecular_mass("BaO"),
        molecular_mass("Li2O"),
        molecular_mass("MgO"),
        molecular_mass("TiO2"),
        molecular_mass("VO2"),
        molecular_mass("Cr2O3"),
        molecular_mass("MnO"),
        molecular_mass("ZnO"),
        molecular_mass("FeO"),
        molecular_mass("Fe2O3"),
        molecular_mass("AlO"),
        molecular_mass("Al2O3"),
        molecular_mass("BeO"),
        molecular_mass("B2O3"),
        molecular_mass("SiO2"),
        atomic_mass('F'),
        atomic_mass("Cl"),
        molecular_mass("H2O")
    ]

    moles_compound = weight_percent ./ molecular_weights
    oxide_factors = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 1, 1, 1, 3, 1, 3, 2, 0.5, 0.5, 0.5]
    moles_oxygen = moles_compound .* oxide_factors
    oxygen_sum = sum(moles_oxygen)
    normalising_factor = normalising_O / oxygen_sum
    moles_oxygen[19] = moles_oxygen[19] * 2
    moles_oxygen[20] = moles_oxygen[20] * 2
    normalised_oxygen_moles = moles_oxygen .* normalising_factor

    # iterate OH
    iterations = 0
    H_initial = 2 - (normalised_oxygen_moles[19] + normalised_oxygen_moles[20])
    println("Initial H: ", H_initial)
    O_OH = (0.5 * H_initial) / normalising_factor
    moles_oxygen[21] = O_OH
    oxygen_sum = sum(moles_oxygen) - 0.5 * (moles_oxygen[19] + moles_oxygen[20])
    normalising_factor = normalising_O / oxygen_sum
    normalised_oxygen_moles = moles_oxygen .* normalising_factor

    ϵ = 1
    while ϵ > 1e-4 && iterations <= 1000
        H_iteration = 2 - (normalised_oxygen_moles[19] + normalised_oxygen_moles[20])
        println("H iteration: ", H_iteration)
        O_OH = (0.5 * H_iteration) / normalising_factor
        moles_oxygen[21] = O_OH
        oxygen_sum = sum(moles_oxygen) - 0.5 * (moles_oxygen[19]+moles_oxygen[20])
        normalising_factor = normalising_O / oxygen_sum
        normalised_oxygen_moles = moles_oxygen .* normalising_factor
        ϵ = abs(H_iteration - H_initial)
        println("Oxygen Sum: ", oxygen_sum, ", Normalising Factor: ", normalising_factor)  # Debugging statement
        iterations = iterations + 1
        H_initial = H_iteration
    end

    element_oxy_ratio = [2, 2, 1, 2, 2, 1, 2, 1, 0.5, 0.5, 0.667, 1, 1, 1, 0.667, 1, 0.667, 0.5, 1, 1, 2]D
    atoms_per_formula_unit = normalised_oxygen_moles .* element_oxy_ratio
    return round.(atoms_per_formula_unit,digits=3)
end

function _find_columns_mica(data)
    workingdata = DataFrame()
    if "sample" in lowercase.(names(data))
        workingdata.Sample = data[:, findfirst(lowercase.(names(data)) .== "sample")]
    else
        workingdata.Sample .= "sample"
    end
    if "Na" in names(data)
        workingdata.Na = data[:, :Na]
    elseif "Na2O" in names(data)
        workingdata.Na = data[:, :Na2O]
    else
        workingdata.Na .= 0
    end
    if "K" in names(data)
        workingdata.K = data[:, :K]
    elseif "K2O" in names(data)
        workingdata.K = data[:, :K2O]
    else
        workingdata.K .= 0
    end
    if "Ca" in names(data)
        workingdata.Ca = data[:, :Ca]
    elseif "CaO" in names(data)
        workingdata.Ca = data[:, :CaO]
    else
        workingdata.Ca .= 0
    end
    if "Rb" in names(data)
        workingdata.Rb = data[:, :Rb]
    elseif "Rb2O" in names(data)
        workingdata.Rb = data[:, :Rb2O]
    else
        workingdata.Rb .= 0
    end
    if "Cs" in names(data)
        workingdata.Cs = data[:, :Cs]
    elseif "Cs2O" in names(data)
        workingdata.Cs = data[:, :Cs2O]
    else
        workingdata.Cs .= 0
    end
    if "Ba" in names(data)
        workingdata.Ba = data[:, :Ba]
    elseif "BaO" in names(data)
        workingdata.Ba = data[:, :BaO]
    else
        workingdata.Ba .= 0
    end
    if "Li" in names(data)
        workingdata.Li = data[:, :Li]
    elseif "Li2O" in names(data)
        workingdata.Li = data[:, :Li2O]
    else
        workingdata.Li .= 0
    end
    if "Mg" in names(data)
        workingdata.Mg = data[:, :Mg]
    elseif "MgO" in names(data)
        workingdata.Mg = data[:, :MgO]
    else
        workingdata.Mg .= 0
    end
    if "Ti" in names(data)
        workingdata.Ti = data[:, :Ti]
    elseif "TiO2" in names(data)
        workingdata.Ti = data[:, :TiO2]
    else
        workingdata.Ti .= 0
    end
    if "V" in names(data)
        workingdata.V = data[:, :V]
    elseif "VO2" in names(data)
        workingdata.V = data[:, :VO2]
    else
        workingdata.V .= 0
    end
    if "Cr" in names(data)
        workingdata.Cr = data[:, :Cr]
    elseif "Cr2O3" in names(data)
        workingdata.Cr = data[:, :Cr2O3]
    else
        workingdata.Cr .= 0
    end
    if "Mn" in names(data)
        workingdata.Mn = data[:, :Mn]
    elseif "MnO" in names(data)
        workingdata.Mn = data[:, :MnO]
    else
        workingdata.Mn .= 0
    end
    if "Zn" in names(data)
        workingdata.Zn = data[:, :Zn]
    elseif "ZnO" in names(data)
        workingdata.Zn = data[:, :ZnO]
    else
        workingdata.Zn .= 0
    end
    if "FeO" in names(data)
            workingdata.FeO = data[:, :FeO]
    else
        workingdata.FeO = 0
    end
    if "Fe2O3" in names(data)
        workingdata.Fe2O3 = data[:, :Fe2O3]
    end
    if "Fe" in names(data)
            workingdata.FeO = 0
            workingdata.Fe2O3 = 0
            workingdata.Fe = data[:, :Fe]
    else
        workingdata.Fe .= 0
    end
    if "Al" in names(data)
        workingdata.Al = data[:, :Al]
    elseif "Al2O3" in names(data)
        workingdata.Al = data[:, :Al2O3]
    else
        workingdata.Al .= 0
    end
    if "Be" in names(data)
        workingdata.Be = data[:, :Be]
    elseif "BeO" in names(data)
        workingdata.Be = data[:, :BeO]
    else
        workingdata.Be .= 0
    end
    if "B" in names(data)
        workingdata.B = data[:, :B]
    elseif "B2O3" in names(data)
        workingdata.B = data[:, :B2O3]
    else
        workingdata.B .= 0
    end
    if "Si" in names(data)
        workingdata.Si = data[:, :Si]
    elseif "SiO2" in names(data)
        workingdata.Si = data[:, :SiO2]
    else
        workingdata.Si .= 0
    end
    if "F" in names(data)
        workingdata.F = data[:, :F]
    else
        workingdata.F .= 0
    end
    if "Cl" in names(data)
        workingdata.Cl = data[:, :Cl]
    else
        workingdata.Cl .= 0
    end
    if "OH" in names(data)
        workingdata.OH = data[:, :OH]
    elseif "H2O" in names(data)
        workingdata.OH = data[:, :H2O]
    else
        workingdata.OH .= 0
    end
    return workingdata
end
