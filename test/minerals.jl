@testset "oxide conversion" begin
    @test element_to_oxide(1; element="Si",element_multiplicity=1,oxide="SiO2") ≈ 2.13933524
    @test element_to_oxide(10000; element="Si",element_multiplicity=1,oxide="SiO2", units="ppm") ≈ 2.13933524
    @test element_to_oxide(1; element="Al",element_multiplicity=2,oxide="Al2O3") ≈ 1.88946373
    @test element_to_oxide(10000; element="Al",element_multiplicity=2,oxide="Al2O3", units="ppm") ≈ 1.88946373
    @test element_to_oxide(1; element="K",element_multiplicity=2,oxide="K2O") ≈ 1.2046048
    @test element_to_oxide(10000; element="K",element_multiplicity=2,oxide="K2O", units="ppm") ≈ 1.2046048
    @test element_to_oxide(1; element='K',element_multiplicity=2,oxide="K2O") ≈ 1.2046048
    @test element_to_oxide(10000; element='K',element_multiplicity=2,oxide="K2O", units="ppm") ≈ 1.2046048
end

@testset "molecular mass" begin
    # @test molecular_mass("Na0.84Ca0.64K0.44Ba0.02Sr0.01Ca0.64Mn0.05Zn0.01Nb2.7Ti1.28Fe0.07Si7.97Al0.06O24.09(OH)1.23O2.75(H2O)6.92") ≈ 1207.6534935612 #Karupmollerite-Ca
    @test molecular_mass(["Al2O3", "K2O"]) ≈ [101.9612772, 94.196]
    @test molecular_mass("(Na3Na3Na3Na3Na3)Ca6(Fe2)Zr3(Si3O9)2(Si9O27SiO)(Si9O27SiO)Cl(H2O)") ≈ 2938.3234192
end
