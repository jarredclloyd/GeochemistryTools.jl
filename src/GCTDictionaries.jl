#=
This jl file contains dictionaries and other reference information for that are used in "GeochemistryTools.jl"
=#

export Avogadro, element_symbol_to_mass, Dict_SrInitial,
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
const Avogadro = 6.02214076 * 10^23

#Constant Declarations: decay constants
const λK40, SEλK40 = 0.00055305, 0.00000132
const λRb87, SEλRb87 = 0.000013972, 0.000000045
const λSm147, SEλSm147 = 0.000006524, 0.000000012
const λLu176, SEλLu176 = 0.00001867, 0.00000008
const λRe187, SEλRe187 = 0.00001666, 0.000000085
const λRa226, SEλRa226 = 0.4332, 0.0019
const λTh230, SEλTh230 = 0.0091705, 0.0000016
const λTh232, SEλTh232 = 0.0000495, 0.0000025
const λPa231, SEλPa231 = 0.021158, 0.00071
const λU234, SEλU234 = 0.00282206, 0.00000080
const λU235, SEλU235 = 0.00098485, 0.000000670
const λU238, SEλU238 = 0.000155125, 0.000000083

#Constant Declarations: isotope ratios
const U238U235, SE238U235 = 137.818, 0.0225

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
element_symbol_to_mass = Dict(
    [('H', 1.00794), ("He", 4.002602), ("Li", 6.941), ("Be", 9.012182), ('B', 10.811), ('C', 12.0107), ('N', 14.0067),
    ('O', 15.9994), ('F', 18.9984032), ("Ne", 20.1797), ("Na", 22.98976928), ("Mg", 24.305), ("Al", 26.9815386),
    ("Si", 28.0855), ('P', 30.973762), ('S', 32.065), ("Cl", 35.453), ("Ar", 39.948), ('K', 39.0983), ("Ca", 40.078),
    ("Sc", 44.955912), ("Ti", 47.867), ('V', 50.9415), ("Cr", 51.9961), ("Mn", 54.938045), ("Fe", 55.845),
    ("Co", 58.933195), ("Ni", 58.6934), ("Cu", 63.546), ("Zn", 65.409), ("Ga", 69.723), ("Ge", 72.64), ("As", 74.9216),
    ("Se", 78.96), ("Br", 79.904), ("Kr", 83.798), ("Rb", 85.4678), ("Sr", 87.62), ('Y', 88.90585), ("Zr", 91.224),
    ("Nb", 92.90638), ("Mo", 95.94), ("Ru", 101.07), ("Rh", 102.9055), ("Pd", 106.42), ("Ag", 107.8682), ("Cd", 112.411),
    ("In", 114.818), ("Sn", 118.71), ("Sb", 121.76), ("Te", 127.6), ('I', 126.90447), ("Xe", 131.293),
    ("Cs", 132.9054519), ("Ba", 137.327), ("La", 138.90547), ("Ce", 140.116), ("Pr", 140.90765), ("Nd", 144.242),
    ("Sm", 150.36), ("Eu", 151.964), ("Gd", 157.25), ("Tb", 158.92535), ("Dy", 162.5), ("Ho", 164.93032),
    ("Er", 167.259), ("Tm", 168.93421), ("Yb", 173.04), ("Lu", 174.967), ("Hf", 178.49), ("Ta", 180.94788),
    ("W", 183.84), ("Re", 186.207), ("Os", 190.23), ("Ir", 192.217), ("Pt", 195.084), ("Au", 196.966569), ("Hg", 200.59),
    ("Tl", 204.3833), ("Pb", 207.2), ("Bi", 208.9804), ("Th", 232.03806), ("Pa", 231.03588), ('U', 238.0289)]
)

#Isotopic properties (mass, radii, charge)

#C1 Chondrite (Palme and O'Neill 2016)

# Initial Sr87/86 ratios
Dict_SrInitial = Dict(
    [
        ("MDC", [1e-10, 1e-12, 0.72607, 0.00070]), ("MDCInv", [1e-10, 1e-12, 1.37727767295, 0.0013278]),
        ("Hogsbo", [1e-10, 1e-12, 0.72000, 0.000103]), ("HogsboInv", [1e-10, 1e-12, 1.38889, 0.000198688]),
        ("MicaMg", [1e-10, 1e-12, 0.72607, 0.00070]), ("MicaMgInv", [1e-10, 1e-12, 1.37727767295, 0.0013278]),
        ("SolarSystem", [1e-10, 1e-12,  0.69899, 0.00005]), ("SolarSystemInv", [1e-10, 1e-12, 1.430635631, 0.000102336]), 
        ("GLO", [1e-10, 1e-12, 0.707358, 0.0001]), ("GLOInv", [1e-10, 1e-12, 1.413727292, 0.0002]),
        ("LaPosta", [1e-10, 1e-12,  0.70483, 0.0001]), ("LaPostaInv", [1e-10, 1e-12, 1.418781834, 0.0002]), 
        ("MtDromedary", [1e-10, 1e-12,  0.7049, 0.0001]), ("MtDromedaryInv", [1e-10, 1e-12, 1.418640942, 0.0002]), 
        ("WilsonsProm", [1e-10, 1e-12, 0.71, 0.0001]), ("WilsonsPromInv", [1e-10, 1e-12, 1.408450704, 0.0002]), 
        ("NIST610", [1e-10, 1e-12, 0.7097, 0.0001]), ("NIST610Inv", [1e-10, 1e-12, 1.409, 0.0002]),
        ("BCR2G", [1e-10, 1e-12, 0.704913, 0.00001]), ("BCR2GInv", [1e-10, 1e-12, 1.418615, 0.000020]),
        ("NIST612", [1e-10, 1e-12, 0.7090630, 0.0001]), ("NIST612Inv", [1e-10, 1e-12, 1.4103120, 0.0002])
    ]
)