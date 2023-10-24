#=
This jl file contains dictionaries and other reference information for that are used in "GeochemistryTools.jl"
=#

export Avogadro,
    element_symbol_to_mass,
    molecular_weight_oxide,
    dict_sr87_sr86i,
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
element_symbol_to_name = Dict()

"""
    element_symbol_to_mass(x)
Return the IUPAC aggregate atomic mass (gmol⁻¹) for a given element.

Single character elements must be entered as a character using ' '.
Double character elements must be entered as a string " ".

#Examples
```julia-repl
julia > element_symbol_to_mass('H')
1.00794
julia> element_symbol_to_mass("He")
4.002602
```
"""
const element_symbol_to_mass::AbstractDict = Dict{Any, AbstractFloat}(
    [('H', 1.00794), ("He", 4.002602), ("Li", 6.941), ("Be", 9.012182), ('B', 10.811),
    ('C', 12.0107), ('N', 14.0067), ('O', 15.9994), ('F', 18.9984032), ("Ne", 20.1797),
    ("Na", 22.98976928), ("Mg", 24.305), ("Al", 26.9815386), ("Si", 28.0855),
    ('P', 30.973762), ('S', 32.065), ("Cl", 35.453), ("Ar", 39.948), ('K', 39.0983),
    ("Ca", 40.078), ("Sc", 44.955912), ("Ti", 47.867), ('V', 50.9415), ("Cr", 51.9961),
    ("Mn", 54.938045), ("Fe", 55.845), ("Co", 58.933195), ("Ni", 58.6934), ("Cu", 63.546),
    ("Zn", 65.409), ("Ga", 69.723), ("Ge", 72.64), ("As", 74.9216), ("Se", 78.96),
    ("Br", 79.904), ("Kr", 83.798), ("Rb", 85.4678), ("Sr", 87.62), ('Y', 88.90585),
    ("Zr", 91.224), ("Nb", 92.90638), ("Mo", 95.94), ("Ru", 101.07), ("Rh", 102.9055),
    ("Pd", 106.42), ("Ag", 107.8682), ("Cd", 112.411), ("In", 114.818), ("Sn", 118.71),
    ("Sb", 121.76), ("Te", 127.6), ('I', 126.90447), ("Xe", 131.293), ("Cs", 132.9054519),
    ("Ba", 137.327), ("La", 138.90547), ("Ce", 140.116), ("Pr", 140.90765), ("Nd", 144.242),
    ("Sm", 150.36), ("Eu", 151.964), ("Gd", 157.25), ("Tb", 158.92535), ("Dy", 162.5),
    ("Ho", 164.93032), ("Er", 167.259), ("Tm", 168.93421), ("Yb", 173.04), ("Lu", 174.967),
    ("Hf", 178.49), ("Ta", 180.94788), ("W", 183.84), ("Re", 186.207), ("Os", 190.23),
    ("Ir", 192.217), ("Pt", 195.084), ("Au", 196.966569), ("Hg", 200.59),
    ("Tl", 204.3833), ("Pb", 207.2), ("Bi", 208.9804), ("Th", 232.03806), ("Pa", 231.03588),
    ('U', 238.0289)]
)

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
]
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

# Initial Sr87/86 ratios
"""
    dict_sr87_sr86i["key"]
Return the inital ⁸⁷Sr/⁸⁶Sr for a given reference material and its standard error

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
