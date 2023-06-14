#= Preamble
This file contains functions for calculating mineral formulas from ppm or wt% data
=#
export MicaFormula

"""
    MicaFormula(data, units; [normalising_O::Int=11])

Calculates mica formula based on input data and declared units.

# Description
Data is a data frame containing sample name and elemental data (either as element or oxide) in a consistent unit
(either ppm OR wt%).

Units is the units the data is in, a string of value "ppm" OR "wt%"

# Example
```julia-repl
julia> MicaFormula(df, "ppm")
```
"""

function MicaFormula(data::AbstractDataFrame, units::AbstractString; normalising_O::Int=11)
    I_site = [:Na, :K, :Ca, :Rb, :Cs, :Ba]
    M_site = [:Li, :Mg, :Ti, :V, :Cr, :Mn, :Zn, :Fe]
    TM_site = [:Al]
    T_site = [:Be, :B, :Si]
    A_site = [:F, :Cl, :OH]
    column_names = [:Na, :K, :Ca, :Rb, :Cs, :Ba, :Li, :Mg, :Ti, :V, :Cr, :Mn, :Zn, :Fe, :Al, :Be, :B, :Si, :F, :Cl, :OH]
    workingdata = _micaFindColumns(data)
    workingvector = collect(values(workingdata[1, Not(:Sample)]))

    if uppercase(units) == "PPM"
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
        molecular_weight_oxide["Na2O"], molecular_weight_oxide["K2O"], molecular_weight_oxide["CaO"],
        molecular_weight_oxide["Rb2O"], molecular_weight_oxide["Cs2O"], molecular_weight_oxide["BaO"],
        molecular_weight_oxide["Li2O"], molecular_weight_oxide["MgO"], molecular_weight_oxide["TiO2"],
        molecular_weight_oxide["VO2"], molecular_weight_oxide["Cr2O3"], molecular_weight_oxide["MnO"],
        molecular_weight_oxide["ZnO"], molecular_weight_oxide["FeO"], molecular_weight_oxide["Al2O3"],
        molecular_weight_oxide["BeO"], molecular_weight_oxide["B2O3"], molecular_weight_oxide["SiO2"],
        element_symbol_to_mass['F'], element_symbol_to_mass["Cl"], molecular_weight_oxide["H2O"]
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
    O_OH = (0.5 * H_initial) / normalising_factor
    moles_oxygen[21] = O_OH
    oxygen_sum = sum(moles_oxygen) - 0.5 * (moles_oxygen[19 + moles_oxygen[20]])
    normalising_factor = normalising_O / oxygen_sum
    normalised_oxygen_moles = moles_oxygen .* normalising_factor

    ϵ = 1
    while normalised_oxygen_sum > normalising_O && ϵ > 1e-8 && iterations <= 1000
        # moles_oxygen = moles_compound .* oxide_factors
        # oxygen_sum = sum(moles_oxygen)
        # normalising_factor = normalising_O / oxygen_sum


        # normalised_oxygen_moles = moles_oxygen .* normalising_factor
        # oxygen_sum = sum(moles_oxygen) - (0.5 * (moles_oxygen[19] + moles_oxygen[20]))
        # normalising_factor = normalising_O / oxygen_sum
        # normalised_oxygen_sum = sum(normalised_oxygen_moles)
        # H_iteration = 2 - (normalised_oxygen_moles[19] + normalised_oxygen_moles[20])
        # O_OH = (0.5 * H_iteration) / normalising_factor
        # moles_oxygen[21] = O_OH
        # normalised_oxygen_moles[21] = O_OH
        # ϵ = abs(H_iteration - H_initial)
        # iterations = iterations + 1
        # H_initial = H_iteration
    end
    element_oxy_ratio = [2, 2, 1, 2, 2, 1, 2, 1, 0.5, 0.5, 0.667, 1, 1, 1, 0.667, 1, 0.667, 0.5, 0.5, 0.5, 0.5]
    atoms_per_formula_unit = normalised_oxygen_moles .* element_oxy_ratio
    return atoms_per_formula_unit, iterations
end

# use lowercase.(names(data)) to make code shorter.
function _micaFindColumns(data)
    workingdata = DataFrame()
    if in("sample", lowercase.(names(data))) == true
        workingdata.Sample = data[:, :Sample]
    else
        workingdata.Sample .= "sample"
    end
    if in("Na", names(data)) == true
        workingdata.Na = data[:, :Na]
    elseif in("Na2O", names(data)) == true
        workingdata.Na = data[:, :Na2O]
    else
        workingdata.Na .= 0
    end
    if in("K", names(data)) == true
        workingdata.K = data[:, :K]
    elseif in("K2O", names(data)) == true
        workingdata.K = data[:, :K2O]
    else
        workingdata.K .= 0
    end
    if in("Ca", names(data)) == true
        workingdata.Ca = data[:, :Ca]
    elseif in("CaO", names(data)) == true
        workingdata.Ca = data[:, :CaO]
    else
        workingdata.Ca .= 0
    end
    if in("Rb", names(data)) == true
        workingdata.Rb = data[:, :Rb]
    elseif in("Rb2O", names(data)) == true
        workingdata.Rb = data[:, :Rb2O]
    else
        workingdata.Rb .= 0
    end
    if in("Cs", names(data)) == true
        workingdata.Cs = data[:, :Cs]
    elseif in("Cs2O", names(data)) == true
        workingdata.Cs = data[:, :Cs2O]
    else
        workingdata.Cs .= 0
    end
    if in("Ba", names(data)) == true
        workingdata.Ba = data[:, :Ba]
    elseif in("BaO", names(data)) == true
        workingdata.Ba = data[:, :BaO]
    else
        workingdata.Ba .= 0
    end
    if in("Li", names(data)) == true
        workingdata.Li = data[:, :Li]
    elseif in("Li2O", names(data)) == true
        workingdata.Li = data[:, :Li2O]
    else
        workingdata.Li .= 0
    end
    if in("Mg", names(data)) == true
        workingdata.Mg = data[:, :Mg]
    elseif in("MgO", names(data)) == true
        workingdata.Mg = data[:, :MgO]
    else
        workingdata.Mg .= 0
    end
    if in("Ti", names(data)) == true
        workingdata.Ti = data[:, :Ti]
    elseif in("TiO2", names(data)) == true
        workingdata.Ti = data[:, :TiO2]
    else
        workingdata.Ti .= 0
    end
    if in("V", names(data)) == true
        workingdata.V = data[:, :V]
    elseif in("VO2", names(data)) == true
        workingdata.V = data[:, :VO2]
    else
        workingdata.V .= 0
    end
    if in("Cr", names(data)) == true
        workingdata.Cr = data[:, :Cr]
    elseif in("Cr2O3", names(data)) == true
        workingdata.Cr = data[:, :Cr2O3]
    else
        workingdata.Cr .= 0
    end
    if in("Mn", names(data)) == true
        workingdata.Mn = data[:, :Mn]
    elseif in("MnO", names(data)) == true
        workingdata.Mn = data[:, :MnO]
    else
        workingdata.Mn .= 0
    end
    if in("Zn", names(data)) == true
        workingdata.Zn = data[:, :Zn]
    elseif in("ZnO", names(data)) == true
        workingdata.Zn = data[:, :ZnO]
    else
        workingdata.Zn .= 0
    end
    if in("Fe", names(data)) == true
        workingdata.Fe = data[:, :Fe]
    elseif in("FeO", names(data)) == true
        workingdata.Fe = data[:, :FeO]
    else
        workingdata.Fe .= 0
    end
    if in("Al", names(data)) == true
        workingdata.Al = data[:, :Al]
    elseif in("Al2O3", names(data)) == true
        workingdata.Al = data[:, :Al2O3]
    else
        workingdata.Al .= 0
    end
    if in("Be", names(data)) == true
        workingdata.Be = data[:, :Be]
    elseif in("BeO", names(data)) == true
        workingdata.Be = data[:, :BeO]
    else
        workingdata.Be .= 0
    end
    if in("B", names(data)) == true
        workingdata.B = data[:, :B]
    elseif in("B2O3", names(data)) == true
        workingdata.B = data[:, :B2O3]
    else
        workingdata.B .= 0
    end
    if in("Si", names(data)) == true
        workingdata.Si = data[:, :Si]
    elseif in("SiO2", names(data)) == true
        workingdata.Si = data[:, :SiO2]
    else
        workingdata.Si .= 0
    end
    if in("F", names(data)) == true
        workingdata.F = data[:, :F]
    else
        workingdata.F .= 0
    end
    if in("Cl", names(data)) == true
        workingdata.Cl = data[:, :Cl]
    else
        workingdata.Cl .= 0
    end
    if in("OH", names(data)) == true
        workingdata.OH = data[:, :OH]
    elseif in("H2O", names(data)) == true
        workingdata.OH = data[:, :H2O]
    else
        workingdata.OH .= 0
    end
    return workingdata
end
