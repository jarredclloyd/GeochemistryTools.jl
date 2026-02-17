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
function formula_mica(
    data::AbstractDataFrame,
    units::AbstractString;
    normalising_O::Int = 11,
    max_iterations::Int = 100,
    tolerance::Float64 = 1e-5,
    Fe_priority::Float64 = 0.7,
)
    I_site = [:NH4, :Na, :K, :Ca, :Rb, :Cs, :Ba]
    M_site = [:Li, :Mg, :Ti, :V, :Cr, :MnO, :Mn2O3, :Zn, :FeO, :Fe2O3, :Al]
    T_site = [:Fe2O3, :Al, :Be, :B, :Si]
    A_site = [:F, :Cl, :OH, :S]
    column_names = [
        :NH4,
        :Na,
        :K,
        :Ca,
        :Rb,
        :Cs,
        :Ba,
        :Li,
        :Mg,
        :Ti,
        :V,
        :Cr,
        :MnO,
        :Mn2O3,
        :Zn,
        :FeO,
        :Fe2O3,
        :Al,
        :Be,
        :B,
        :Si,
        :F,
        :Cl,
        :OH,
        :S,
    ]
    workingdata = _find_columns_mica(data)
    workingvector = collect(values(workingdata[1, Not(:Sample)]))

    if occursin("ppm", lowercase(units))
        oxide_conversion_factors = [
            cation_to_oxide(1; cation = "NH4", cation_multiplicity = 2, oxide = "(NH4)2O"),
            element_to_oxide(
                1;
                element = :Na,
                element_multiplicity = 2,
                oxide = "Na2O",
            ).val,
            element_to_oxide(1; element = :K, element_multiplicity = 2, oxide = "K2O").val,
            element_to_oxide(1; element = :Ca, element_multiplicity = 1, oxide = "CaO").val,
            element_to_oxide(
                1;
                element = :Rb,
                element_multiplicity = 2,
                oxide = "Rb2O",
            ).val,
            element_to_oxide(
                1;
                element = :Cs,
                element_multiplicity = 2,
                oxide = "Cs2O",
            ).val,
            element_to_oxide(1; element = :Ba, element_multiplicity = 1, oxide = "BaO").val,
            element_to_oxide(
                1;
                element = :Li,
                element_multiplicity = 2,
                oxide = "Li2O",
            ).val,
            element_to_oxide(1; element = :Mg, element_multiplicity = 1, oxide = "MgO").val,
            element_to_oxide(
                1;
                element = :Ti,
                element_multiplicity = 1,
                oxide = "TiO2",
            ).val,
            element_to_oxide(1; element = :V, element_multiplicity = 2, oxide = "V2O3").val,
            element_to_oxide(
                1;
                element = :Cr,
                element_multiplicity = 2,
                oxide = "Cr2O3",
            ).val,
            element_to_oxide(1; element = :Mn, element_multiplicity = 1, oxide = "MnO").val,
            element_to_oxide(
                1;
                element = :Mn,
                element_multiplicity = 2,
                oxide = "Mn2O3",
            ).val,
            element_to_oxide(1; element = :Zn, element_multiplicity = 1, oxide = "ZnO").val,
            element_to_oxide(1; element = :Fe, element_multiplicity = 1, oxide = "FeO").val,
            element_to_oxide(
                1;
                element = :Fe,
                element_multiplicity = 2,
                oxide = "Fe2O3",
            ).val,
            element_to_oxide(
                1;
                element = :Al,
                element_multiplicity = 2,
                oxide = "Al2O3",
            ).val,
            element_to_oxide(1; element = :Be, element_multiplicity = 1, oxide = "BeO").val,
            element_to_oxide(1; element = :B, element_multiplicity = 2, oxide = "B2O3").val,
            element_to_oxide(
                1;
                element = :Si,
                element_multiplicity = 1,
                oxide = "SiO2",
            ).val,
            element_to_oxide(1; element = :F, element_multiplicity = 1, oxide = "F").val,
            element_to_oxide(1; element = :Cl, element_multiplicity = 1, oxide = "Cl").val,
            1, #OH
            element_to_oxide(1; element = :S, element_multiplicity = 1, oxide = "S").val,
        ]
        weight_percent = (workingvector ./ 1e4) .* oxide_conversion_factors
    elseif occursin("wt", lowercase(units))
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
        molecular_mass("V2O3"; verbose = false),
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
        get_atomicmass(:S).val,
    ]

    moles_compound = weight_percent ./ molecular_weights
    oxy_per_elem = [
        0.5,# NH4 (as (NH4)₂O equivalent)
        1.0,  # Na (as Na₂O)
        1.0,  # K (as K₂O)
        1.0,  # Ca (as CaO)
        1.0,  # Rb (as Rb₂O)
        1.0,  # Cs (as Cs₂O)
        1.0,  # Ba (as BaO)
        1.0,  # Li (as Li₂O)
        1.0,  # Mg (as MgO)
        2.0,  # Ti (as TiO₂)
        2.5,  # V (as V₂O₅)
        1.5,  # Cr (as Cr₂O₃)
        1.0,  # MnO (as MnO)
        1.5,  # Mn₂O₃ (as Mn₂O₃)
        1.0,  # Zn (as ZnO)
        1.0,  # FeO (as FeO)
        1.5,  # Fe₂O₃ (as Fe₂O₃)
        1.5,  # Al (as Al₂O₃)
        1.0,  # Be (as BeO)
        1.5,  # B (as B₂O₃)
        2.0,  # Si (as SiO₂)
        0.0,  # F (substitutes for OH)
        0.0,  # Cl (substitutes for OH)
        1.0,  # OH
        0.0,   # S (as S²⁻, substitutes for O²⁻)
    ]
    cation_per_elem = [
        2.0,  # NH4 (as (NH4)₂O)
        2.0,  # Na (as Na₂O)
        2.0,  # K (as K₂O)
        1.0,  # Ca (as CaO)
        2.0,  # Rb (as Rb₂O)
        2.0,  # Cs (as Cs₂O)
        1.0,  # Ba (as BaO)
        2.0,  # Li (as Li₂O)
        1.0,  # Mg (as MgO)
        1.0,  # Ti (as TiO₂)
        2.0,  # V (as V₂O₅)
        2.0,  # Cr (as Cr₂O₃)
        1.0,  # MnO (as MnO)
        2.0,  # Mn₂O₃ (as Mn₂O₃)
        1.0,  # Zn (as ZnO)
        1.0,  # FeO (as FeO)
        2.0,  # Fe₂O₃ (as Fe₂O₃)
        2.0,  # Al (as Al₂O₃)
        1.0,  # Be (as BeO)
        2.0,  # B (as B₂O₃)
        1.0,  # Si (as SiO₂)
        1.0,  # F
        1.0,  # Cl
        1.0,  # OH
        1.0,   # S
    ]
    oxidation_states = [
        1.0,  # NH4⁺
        1.0,  # Na⁺
        1.0,  # K⁺
        2.0,  # Ca²⁺
        1.0,  # Rb⁺
        1.0,  # Cs⁺
        2.0,  # Ba²⁺
        1.0,  # Li⁺
        2.0,  # Mg²⁺
        4.0,  # Ti⁴⁺
        5.0,  # V⁵⁺
        3.0,  # Cr³⁺
        2.0,  # Mn²⁺ (MnO)
        3.0,  # Mn³⁺ (Mn₂O₃)
        2.0,  # Zn²⁺
        2.0,  # Fe²⁺ (FeO)
        3.0,  # Fe³⁺ (Fe₂O₃)
        3.0,  # Al³⁺
        2.0,  # Be²⁺
        3.0,  # B³⁺
        4.0,  # Si⁴⁺
        -1.0, # F⁻
        -1.0, # Cl⁻
        -1.0, # OH⁻
        -2.0,  # S²⁻
    ]
    moles_oxygen = moles_compound .* oxy_per_elem
    oxygen_sum = sum(moles_oxygen)
    normalising_factor = normalising_O / oxygen_sum
    normalised_moles = moles_compound .* normalising_factor

    cation_formula = normalised_moles .* cation_per_elem

    # normalised_oxygen_moles = moles_oxygen .* normalising_factor

    # iterate OH
    OH = cation_formula[24]
    ϵ = Inf

    if OH == 0
        iteration = 0
        OH_est = 2 - (cation_formula[22] + cation_formula[23])
        OH_est = clamp(OH_est, 0, 2)
        while ϵ > tolerance && iteration ≤ max_iterations
            iteration += 1
            OH_prior = OH_est
            moles_oxygen = cation_formula ./ cation_per_elem .* oxy_per_elem
            oxygen_sum = sum(moles_oxygen)
            normalising_factor = normalising_O / oxygen_sum
            cation_formula .*= normalising_factor

            OH_est = 2 - (cation_formula[22] + cation_formula[23])
            OH_est = clamp(OH_est, 0, 2)
            println("H iteration: $OH_est")
            cation_formula[24] = OH_est

            # Charge balance
            positive_charge = sum(cation_formula[1:21] .* oxidation_states[1:21])
            required_charge = 22 - (cation_formula[22] + cation_formula[23])
            charge_imbalance = required_charge - positive_charge

            total_Mn = cation_formula[13] + 2 * cation_formula[14]
            total_Fe = cation_formula[16] + 2 * cation_formula[17]
            Fe3 = clamp(charge_imbalance * Fe_priority, 0.0, total_Fe)
            Mn3 = clamp(charge_imbalance * (1.0 - Fe_priority), 0.0, total_Mn)
            cation_formula[13] = total_Mn - Mn3
            cation_formula[14] = Mn3 / 2.0
            cation_formula[16] = total_Fe - Fe3
            cation_formula[17] = Fe3 / 2.0

            ϵ = abs((OH_est - OH_prior) / OH_prior)
            println("epsilon: $ϵ")
        end

        moles_oxygen = cation_formula ./ cation_per_elem .* oxy_per_elem
        oxygen_sum = sum(moles_oxygen)
        normalising_factor = normalising_O / oxygen_sum
        cation_formula .*= normalising_factor
    end

    positive_charge = sum(cation_formula[1:21] .* oxidation_states[1:21])
    charge_imbalance =
        abs(positive_charge - (22 - (cation_formula[22] + cation_formula[23])))


    # Site allocation
    I_site = sum(cation_formula[1:7])
    T_site = sum(cation_formula[19:21])
    M_site = sum(cation_formula[8:16])
    A_site = sum(cation_formula[22:end])
    if T_site < 4.0
        Al = cation_formula[18]
        Al_T = clamp(4 - T_site, 0, 4)
        Al_M = clamp(Al - Al_T, 0, Inf64)
        T_site += Al_T
        M_site += Al_M
        if T_site < 4.0
            Fe3 = cation_formula[17]
            Fe3_T = clamp(4 - T_site, 0, 4)
            Fe3_M = clamp(Fe3 - Fe3_T, 0, Inf)
            T_site += Fe3_T
            M_site += Fe3_M
        end
    end

    sites = Dict(:I => I_site, :T => T_site, :M => M_site, :A => A_site)

    println("Cation sum: $(sum(cation_formula[1:21]))")
    println("Total ions: $(sum(cation_formula))")
    return (
        column_names,
        round.(cation_formula; digits = 3),
        charge_imbalance,
        iteration,
        sites,
    )
end

function _find_columns_mica(data)
    workingdata = DataFrame()
    if in("sample", lowercase.(names(data)))
        workingdata.Sample = data[:, findfirst(lowercase.(names(data)) .== "sample")]
    else
        workingdata.Sample = collect(1:1:size(data, 1))
    end
    if in("NH4", names(data))
        workingdata.NH4 = data[:, :NH4]
    else
        workingdata.NH4 .= 0
    end
    if in("Na", names(data))
        workingdata.Na = data[:, :Na]
    elseif in("Na2O", names(data))
        workingdata.Na = data[:, :Na2O]
    else
        workingdata.Na .= 0
    end
    if in("K", names(data))
        workingdata.K = data[:, :K]
    elseif in("K2O", names(data))
        workingdata.K = data[:, :K2O]
    else
        workingdata.K .= 0
    end
    if in("Ca", names(data))
        workingdata.Ca = data[:, :Ca]
    elseif in("CaO", names(data))
        workingdata.Ca = data[:, :CaO]
    else
        workingdata.Ca .= 0
    end
    if in("Rb", names(data))
        workingdata.Rb = data[:, :Rb]
    elseif in("Rb2O", names(data))
        workingdata.Rb = data[:, :Rb2O]
    else
        workingdata.Rb .= 0
    end
    if in("Cs", names(data))
        workingdata.Cs = data[:, :Cs]
    elseif in("Cs2O", names(data))
        workingdata.Cs = data[:, :Cs2O]
    else
        workingdata.Cs .= 0
    end
    if in("Ba", names(data))
        workingdata.Ba = data[:, :Ba]
    elseif in("BaO", names(data))
        workingdata.Ba = data[:, :BaO]
    else
        workingdata.Ba .= 0
    end
    if in("Li", names(data))
        workingdata.Li = data[:, :Li]
    elseif in("Li2O", names(data))
        workingdata.Li = data[:, :Li2O]
    else
        workingdata.Li .= 0
    end
    if in("Mg", names(data))
        workingdata.Mg = data[:, :Mg]
    elseif in("MgO", names(data))
        workingdata.Mg = data[:, :MgO]
    else
        workingdata.Mg .= 0
    end
    if in("Ti", names(data))
        workingdata.Ti = data[:, :Ti]
    elseif in("TiO2", names(data))
        workingdata.Ti = data[:, :TiO2]
    else
        workingdata.Ti .= 0
    end
    if in("V", names(data))
        workingdata.V = data[:, :V]
    elseif in("VO2", names(data))
        workingdata.V = data[:, :VO2]
    else
        workingdata.V .= 0
    end
    if in("Cr", names(data))
        workingdata.Cr = data[:, :Cr]
    elseif in("Cr2O3", names(data))
        workingdata.Cr = data[:, :Cr2O3]
    else
        workingdata.Cr .= 0
    end
    if in("MnO", names(data))
        workingdata.MnO = data[:, :MnO]
    else
        workingdata.MnO .= 0
    end
    if in("Mn2O3", names(data))
        workingdata.Mn2O3 = data[:, :Mn2O3]
    else
        workingdata.Mn2O3 .= 0
    end
    if in("Mn", names(data))
        workingdata.MnO .= 0
    end
    if in("Zn", names(data))
        workingdata.Zn = data[:, :Zn]
    elseif in("ZnO", names(data))
        workingdata.Zn = data[:, :ZnO]
    else
        workingdata.Zn .= 0
    end
    if in("FeO", names(data))
        workingdata.FeO = data[:, :FeO]
    else
        workingdata.FeO .= 0
    end
    if in("Fe2O3", names(data))
        workingdata.Fe2O3 = data[:, :Fe2O3]
    else
        workingdata.Fe2O3 .= 0
    end
    if in("Fe", names(data))
        workingdata.FeO = data[:, :Fe]
    end
    if in("Al", names(data))
        workingdata.Al = data[:, :Al]
    elseif in("Al2O3", names(data))
        workingdata.Al = data[:, :Al2O3]
    else
        workingdata.Al .= 0
    end
    if in("Be", names(data))
        workingdata.Be = data[:, :Be]
    elseif in("BeO", names(data))
        workingdata.Be = data[:, :BeO]
    else
        workingdata.Be .= 0
    end
    if in("B", names(data))
        workingdata.B = data[:, :B]
    elseif in("B2O3", names(data))
        workingdata.B = data[:, :B2O3]
    else
        workingdata.B .= 0
    end
    if in("Si", names(data))
        workingdata.Si = data[:, :Si]
    elseif in("SiO2", names(data))
        workingdata.Si = data[:, :SiO2]
    else
        workingdata.Si .= 0
    end
    if in("F", names(data))
        workingdata.F = data[:, :F]
    else
        workingdata.F .= 0
    end
    if in("Cl", names(data))
        workingdata.Cl = data[:, :Cl]
    else
        workingdata.Cl .= 0
    end
    if in("OH", names(data))
        workingdata.OH = data[:, :OH]
    else
        workingdata.OH .= 0
    end
    if in("S", names(data))
        workingdata.S = data[:, :S]
    else
        workingdata.S .= 0
    end
    return workingdata
end

SSP18 = DataFrame([
    :SiO2 => 34.79,
    :TiO2 => 3.26,
    :Al2O3 => 18.82,
    :FeO => 21.39,
    :MnO => 0.51,
    :MgO => 7.62,
    :CaO => 0,
    :Na2O => 0.12,
    :K2O => 9.66,
    :BaO => 0.14,
    :F => 0.17,
    :Cl => 0.05,
])
