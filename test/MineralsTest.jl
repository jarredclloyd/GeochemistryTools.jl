@testset "oxide conversion" begin
    @test element_to_oxide(1; element="Si",element_multiplicity=1,oxide="SiO2").val ≈ 2.13932704
    @test element_to_oxide(10000; element="Si",element_multiplicity=1,oxide="SiO2", units="ppm").val ≈ 2.13932704
    @test element_to_oxide(1; element="Al",element_multiplicity=2,oxide="Al2O3").val ≈ 1.889441500
    @test element_to_oxide(10000; element="Al",element_multiplicity=2,oxide="Al2O3", units="ppm").val ≈ 1.889441500
    @test element_to_oxide(1; element="K",element_multiplicity=2,oxide="K2O").val ≈ 1.20459968
    @test element_to_oxide(10000; element="K",element_multiplicity=2,oxide="K2O", units="ppm").val ≈ 1.20459968
    @test element_to_oxide(1; element='K',element_multiplicity=2,oxide="K2O").val ≈ 1.20459968
    @test element_to_oxide(10000; element='K',element_multiplicity=2,oxide="K2O", units="ppm").val ≈ 1.20459968
end

@testset "molecular mass" begin
    @test molecular_mass("Na0.84Ca0.64K0.44Ba0.02Sr0.01Ca0.64Mn0.05Zn0.01Nb2.7Ti1.28Fe0.07Si7.97Al0.06O24.09(OH)1.23O2.75(H2O)6.92"; verbose = false) ≈ 1211.31596 #Karupmollerite-Ca
    @test molecular_mass(["Al2O3", "K2O"]; verbose = false) ≈ [101.9600768, 94.1956]
    @test molecular_mass("(Na3Na3Na3Na3Na3)Ca6(Fe2)Zr3(Si3O9)2(Si9O27SiO)(Si9O27SiO)Cl(H2O)"; verbose = false) ≈ 2938.2775392 #Eudialyte
    @test molecular_mass(["CaFe2O4",
            "CaFe2(SiO4)2(PO4)2",
            "Ca6Be4Mn(SiO4)2(Si2O7)2(OH)2",
            "KNa(MoO2)(SO4)2",
            "Ba(CrO4)",
            "NaCa2(Fe4Fe)(Si6Al2)O22(OH)2",
            "AgTlPbAs2S5",
            "NaNaCa(CuFe)(AsO4)3",
            "Ca3SiO5",
            "Ni9BiSbS8",
            "Fe3Mg24Zn18(SO4)4(CO3)2(OH)81",
            "MnS2",
            "MnMn2O4",
            "Na3Ca(Si3Al3)O12(SO4)",
            "CdS",
            "[Pb(H2O)10][Zn12(OH)20(H2O)(SO4)3]",
            "Cu5Ag11Pb76Sb71As17(As)8S224",
            "Ba[Ti3Cr4Fe2Fe2Mg]O19",
            "K(H2O)Mn2(Fe2Ti)(PO4)4[O(OH)](H2O)10(H2O)4",
            "K0.6Cu18[AsO2(OH)2]4[AsO3OH]10(AsO4)(OH)9.6(H2O)18.6"
        ]; verbose = false) ≈ [215.764,
            525.869523996,
            885.9567754,
            382.14806928,
            253.3191,
            990.83584608,
            829.59139,
            622.20132356,
            228.314,
            1115.461,
            3809.502,
            119.058043,
            228.810129,
            562.290923039,
            144.474,
            1818.233,
            34950.9200749,
            1040.578399,
            991.641433992,
            3767.541105
        ]
end

SSP18 = DataFrame(
    [
    :SiO2=>34.79,
    :TiO2=>3.26,
    :Al2O3=>18.82,
    :FeO=>21.39,
    :MnO=>0.51,
    :MgO=>7.62,
    :CaO=>0,
    :Na2O=>0.12,
    :K2O=>9.66,
    :BaO=>0.14,
    :F=>0.17,
    :Cl=>0.05
    ]
    )
G17560 = DataFrame(
        [
        "sample" => "G17560",
        "Na" => 2513.6,
        "K" => 72284.9,
        "Ca" => 204.3,
        "Rb" => 14176.2,
        "Cs" => 2349.4,
        "Ba" => 2.7,
        "Li" => 21683.4,
        "Mg" => 0.8,
        "Ti" => 45.7,
        "V" => 0.0,
        "Cr" => 0.1,
        "Mn" => 1960.8,
        "Zn" => 23.7,
        "Fe" => 6.1,
        "Al" => 143414.0,
        "Be" => 27.2,
        "B" => 126.7,
        "Si" => 228274.4,
        "H" => 0,
        "F" => 65472.5,
        "Cl" => 62.1,
        "O" => 407300.1,
    ]
    )

@testset "Mineral Formula Calculations" begin
    @test formula_mica(SSP18, "wt%O")
end
