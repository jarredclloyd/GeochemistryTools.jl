using GeochemistryTools
using Test

@testset "GeochemistryTools.jl" begin
    molecularwt("Na0.84Ca0.64K0.44Ba0.02Sr0.01Ca0.64Mn0.05Zn0.01Nb2.7Ti1.28Fe0.07Si7.97Al0.06O24.09OH1.23O2.75(H2O)6.92")
    molecularwt(["Al", "MnS", "PbTe", "Ni3Fe", "As8S9", "Na2Ca(UO2)(CO3)3(H2O)6", "Ag6(Ag4Zn2)Sb4S12S", "Na0.84Ca0.64K0.44Ba0.02Sr0.01Ca0.64Mn0.05Zn0.01Nb2.7Ti1.28Fe0.07Si7.97Al0.06O24.09OH1.23O2.75(H2O)6.92"])
end
