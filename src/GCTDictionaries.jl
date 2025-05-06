#=
This jl file contains dictionaries and other reference information for that are used in "GeochemistryTools.jl"
=#

export avogadro,
    molecular_weight_oxide,
    dict_sr87_sr86i, cn_eight_IR,
    ci_chondrite_PO2016,
    #decay constants:
    λK40, SEλK40,
    λRb87, SEλRb87,
    λSm147, SEλSm147,
    λLu176, SEλLu176,
    λRe187, SEλRe187,
    λRa226, SEλRa226,
    λTh230, SEλTh230,
    λTh232, SEλTh232,
    λPa231, SEλPa231,
    λU234, SEλU234,
    λU235, SEλU235,
    λU238, SEλU238,
    #Isotope ratios:
    U238U235, SE238U235

#Constant Declarations: general
const avogadro::AbstractFloat = 6.02214076 * 10^23
const planck::AbstractFloat = 6.62607015e-34 #JHz⁻¹
const speedoflight::AbstractFloat = 299792458 #ms⁻¹
const boltzmann::AbstractFloat = 1.380649e-23 #JK⁻¹

#Constant Declarations: decay constants
const λK40::AbstractFloat, SEλK40::AbstractFloat = 0.00055305, 0.00000132
const λRb87::AbstractFloat, SEλRb87::AbstractFloat = 0.000013972, 0.000000045
const λSm147::AbstractFloat, SEλSm147::AbstractFloat = 0.000006524, 0.000000012
const λLu176::AbstractFloat, SEλLu176::AbstractFloat = 0.00001867, 0.00000008
const λRe187::AbstractFloat, SEλRe187::AbstractFloat = 0.00001666, 0.000000085
const λRa226::AbstractFloat, SEλRa226::AbstractFloat = 0.4332, 0.0019
const λTh230::AbstractFloat, SEλTh230::AbstractFloat = 0.0091705, 0.0000016
const λTh232::AbstractFloat, SEλTh232::AbstractFloat = 0.0000495, 0.0000025
const λPa231::AbstractFloat, SEλPa231::AbstractFloat = 0.021158, 0.00071
const λU234::AbstractFloat, SEλU234::AbstractFloat = 0.00282206, 0.00000080
const λU235::AbstractFloat, SEλU235::AbstractFloat = 0.00098485, 0.000000670
const λU238::AbstractFloat, SEλU238::AbstractFloat = 0.000155125, 0.000000083

#Constant Declarations: isotope ratios
const U238U235::AbstractFloat, SE238U235::AbstractFloat = 137.818, 0.0225

#Element properties (mass, ionic radii, charges, symbol, name)



"""
    cn_eight_IR(x)
Return the eightfold coordinate ionic radius of specified element.

Single character elements must be entered as a character using ' '.
Double character elements must be entered as a string " ".

#Examples
```julia-repl
julia >
```
"""
cn_eight_IR = Dict(
    [
("La_3+", 1.160),
("Ce_3+", 1.143),
("Pr_3+", 1.126),
("Nd_3+", 1.109),
("Pm_3+", 1.093),
("Sm_3+", 1.079),
("Eu_3+", 1.066),
("Gd_3+", 1.053),
("Tb_3+", 1.040),
("Dy_3+", 1.027),
("Ho_3+", 1.015),
("Er_3+", 1.004),
("Tm_3+", 0.994),
("Yb_3+", 0.985),
("Lu_3+", 0.977)]
)

# Molecular properties
const molecular_weight_oxide::AbstractDict = Dict{AbstractString, AbstractFloat}(
    [
    ("H2O", 18.015), ("Li2O" , 29.879), ("BeO" , 25.0112), ("B2O3" , 69.617), ("Na2O" , 61.979), ("MgO" , 40.304),
    ("Al2O3" , 101.961), ("SiO2" , 60.083), ("P2O5" , 141.943), ("SO2" , 64.058), ("SO3" , 80.057),
    ("SO4" , 96.056), ("K2O" , 94.195), ("CaO" , 56.077), ("TiO2" , 79.865), ("V2O5" , 181.879),
    ("VO2" , 82.94), ("Cr2O3" , 151.989), ("Mn2O3" , 157.873), ("MnO" , 70.937), ("MnO2" , 86.936),
    ("Fe2O3" , 159.687), ("Fe3O4" , 231.531), ("FeO" , 71.844), ("Co3O4" , 240.795), ("CoO" , 74.932),
    ("NiO" , 74.692), ("Cu2O" , 143.091), ("CuO" , 79.545), ("ZnO" , 81.379), ("Ga2O3" , 187.443),
    ("GeO2" , 104.628), ("As2O3" , 197.841), ("SeO2" , 110.969), ("Rb2O" , 186.935), ("SrO" , 103.619),
    ("Y2O3" , 225.809), ("ZrO2" , 123.222), ("Nb2O5" , 265.807), ("MoO3" , 143.947), ("Ag2O" , 231.739),
    ("CdO" , 128.409), ("In2O3" , 277.637), ("SnO2" , 150.708), ("Sb2O3" , 291.517), ("TeO2" , 159.598),
    ("Cs2O" , 281.819), ("BaO" , 153.329), ("La2O3" , 325.817), ("Ce2O3" , 328.237), ("CeO2" , 172.118),
    ("Pr2O3" , 329.817), ("Nd2O3" , 336.477), ("Sm2O3" , 348.717), ("Gd2O3" , 362.497), ("Tb2O3" , 365.857), ("Er2O3" , 382.517), ("HfO2" , 210.488), ("Ta2O5" , 441.895), ("WO3" , 231.837), ("HgO" , 216.589),
    ("Tl2O3" , 456.757), ("PbO" , 223.199), ("PbO2" , 239.198), ("Bi2O3" , 465.957), ("ThO2" , 264.038),
    ("U3O8" , 842.082), ("UO2" , 270.028)
]
)
#Isotopic properties (mass, radii, charge)

#C1 Chondrite (Palme and O'Neill 2016)
ci_chondrite_PO2016::AbstractDict = Dict([
('H', [19700,1970]),
("Li", [1.45,0.145]),
("Be", [0.0219,0.001533]),
('B', [0.775,0.0775]),
('C', [34800,3480]),
('N', [2950,442.5]),
('O', [459000,45900]),
('F', [58.2,9.312]),
("Na", [4962,446.58]),
("Mg", [95400,3816]),
("Al", [8400,504]),
("Si", [107000,3210]),
('P', [985,78.8]),
('S', [53500,2675]),
("Cl", [698,104.7]),
('K', [546,49.14]),
("Ca", [9110,546.6]),
("Sc", [5.81,0.3486]),
("Ti", [447,31.29]),
('V', [54.6,3.276]),
("Cr", [2623,131.15]),
("Mn", [1916,114.96]),
("Fe", [186600,7464]),
("Co", [513,20.52]),
("Ni", [10910,763.7]),
("Cu", [133,18.62]),
("Zn", [309,12.36]),
("Ga", [9.62,0.5772]),
("Ge", [32.6,2.934]),
("As", [1.74,0.1566]),
("Se", [20.3,1.421]),
("Br", [3.26,0.489]),
("Rb", [2.32,0.1856]),
("Sr", [7.79,0.5453]),
('Y', [1.46,0.073]),
("Zr", [3.63,0.1815]),
("Nb", [0.283,0.0283]),
("Mo", [0.961,0.0961]),
("Ru", [0.69,0.0345]),
("Rh", [0.132,0.0066]),
("Pd", [0.56,0.0224]),
("Ag", [0.201,0.01809]),
("Cd", [0.674,0.04718]),
("In", [0.0778,0.00389]),
("Sn", [1.63,0.2445]),
("Sb", [0.145,0.0203]),
("Te", [2.28,0.1596]),
('I', [0.53,0.106]),
("Cs", [0.188,0.01128]),
("Ba", [2.42,0.121]),
("La", [0.2414,0.007242]),
("Ce", [0.6194,0.018582]),
("Pr", [0.0939,0.002817]),
("Nd", [0.4737,0.014211]),
("Sm", [0.1536,0.004608]),
("Eu", [0.05883,0.0017649]),
("Gd", [0.2069,0.006207]),
("Tb", [0.03797,0.0011391]),
("Dy", [0.2558,0.007674]),
("Ho", [0.05644,0.0016932]),
("Er", [0.1655,0.004965]),
("Tm", [0.02609,0.0007827]),
("Yb", [0.1687,0.005061]),
("Lu", [0.02503,0.0007509]),
("Hf", [0.1065,0.003195]),
("Ta", [0.015,0.0015]),
('W', [0.096,0.0096]),
("Re", [0.04,0.002]),
("Os", [0.495,0.02475]),
("Ir", [0.469,0.02345]),
("Pt", [0.925,0.04625]),
("Au", [0.148,0.01776]),
("Hg", [0.35,0.175]),
("Tl", [0.14,0.0154]),
("Pb", [2.62,0.2096]),
("Bi", [0.11,0.0099]),
("Th", [0.03,0.0021]),
('U', [0.0081,0.000567]),
]
)

# Initial Sr87/86 ratios
"""
    dict_sr87_sr86i["key"]
Return the initial ⁸⁷Sr/⁸⁶Sr for a given reference material and its standard error

type 'keys(dict_sr87_sr86i)' for a list of available materials
"""
const dict_sr87_sr86i::AbstractDict = Dict{AbstractString, AbstractArray}(
    [
    ("MDC", [0.72607, 0.00035]), ("MDCInv", [1.37727767295164, 0.00066391282594388300]),
    ("Hogsbo", [0.72000, 0.0000515]), ("HogsboInv", [1.38888888888888, 0.00009934413580246910]),
    ("MicaMg", [0.72607, 0.00035]), ("MicaMgInv", [1.37727767295164, 0.00066391282594388300]),
    ("SolarSystem", [0.69899, 0.00005]), ("SolarSystemInv", [1.43063563141103, 0.00010233591549314200]),
    ("GLO", [0.707358, 0.00005]), ("GLOInv", [1.41371130318735, 0.00009992898243798420]),
    ("LaPosta", [0.70483, 0.00005]), ("LaPostaInv", [1.41878183391739, 0.00010064709461270000]),
    ("MtDromedary", [0.7049, 0.00005]), ("MtDromedaryInv", [1.41864094197758, 0.00010062710611275200]),
    ("WilsonsProm", [0.71, 0.00005]), ("WilsonsPromInv", [1.40845070422535, 0.00009918666931164450]),
    ("NIST610", [0.709699, 0.000009]), ("NIST610Inv", [1.40904806122032, 0.00001786874794945860]),
    ("BCR2G", [0.705003, 0.000004]), ("BCR2GInv", [1.41843368042405, 0.00000804781642304532]),
    ("NIST612", [0.7090630, 0.00001]), ("NIST612Inv", [1.41031191868705, 0.00001988979707990760])
]
)
