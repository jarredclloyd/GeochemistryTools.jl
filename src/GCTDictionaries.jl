#=
This jl file contains dictionaries and other reference information for that are used in "GeochemistryTools.jl"
=#

export AVOGADRO,
    INITIAL_Sr,
    CI_CHONDRITE_PO2016,
    IONIC_RADIUS_CNEIGHT,
    #decay constants:
    LAMBDA_K40, LAMBDA_K40_SE,
    LAMBDA_Rb87, LAMBDA_Rb87_SE,
    LAMBDA_Sm147, LAMBDA_Sm147_SE,
    LAMBDA_Lu176, LAMBDA_Lu176_SE,
    LAMBDA_Re187, LAMBDA_Re187_SE,
    LAMBDA_Ra226, LAMBDA_Ra226_SE,
    LAMBDA_Th230, LAMBDA_Th230_SE,
    LAMBDA_Th232, LAMBDA_Th232_SE,
    LAMBDA_Pa231, LAMBDA_Pa231_SE,
    LAMBDA_U234, LAMBDA_U234_SE,
    LAMBDA_U235, LAMBDA_U235_SE,
    LAMBDA_U238, LAMBDA_U238_SE,
    #Isotope ratios:
    RATIO_U238U235, RATIO_U238U235SE

#Constant Declarations: general
const AVOGADRO::Float64x2 = 6.02214076 * 10^23
const PLANCK::Float64x2 = 6.62607015e-34 #JHz⁻¹
const LIGHTSPEED::Integer = 299792458 #ms⁻¹
const BOLTZMANN::Float64x2 = 1.380649e-23 #JK⁻¹

#Constant Declarations: decay constants
const LAMBDA_K40::Float64, LAMBDA_K40_SE::Float64 = 0.00055305, 0.00000132
const LAMBDA_Rb87::Float64, LAMBDA_Rb87_SE::Float64 = 0.000013972, 0.000000045
const LAMBDA_Sm147::Float64, LAMBDA_Sm147_SE::Float64 = 0.000006524, 0.000000012
const LAMBDA_Lu176::Float64, LAMBDA_Lu176_SE::Float64 = 0.00001867, 0.00000008
const LAMBDA_Re187::Float64, LAMBDA_Re187_SE::Float64 = 0.00001666, 0.000000085
const LAMBDA_Ra226::Float64, LAMBDA_Ra226_SE::Float64 = 0.4332, 0.0019
const LAMBDA_Th230::Float64, LAMBDA_Th230_SE::Float64 = 0.0091705, 0.0000016
const LAMBDA_Th232::Float64, LAMBDA_Th232_SE::Float64 = 0.0000495, 0.0000025
const LAMBDA_Pa231::Float64, LAMBDA_Pa231_SE::Float64 = 0.021158, 0.00071
const LAMBDA_U234::Float64, LAMBDA_U234_SE::Float64 = 0.00282206, 0.00000080
const LAMBDA_U235::Float64, LAMBDA_U235_SE::Float64 = 0.00098485, 0.000000670
const LAMBDA_U238::Float64, LAMBDA_U238_SE::Float64 = 0.000155125, 0.000000083

#Constant Declarations: isotope ratios
const RATIO_U238U235::Float64, RATIO_U238U235SE::Float64 = 137.818, 0.0225

#Element properties (mass, ionic radii, charges, symbol, name)

"""
    IONIC_RADIUS_CNEIGHT(x)
Return the eightfold coordinate ionic radius of specified element.

Request value by key as string with the element symbol followed by an underscore and charge state e.g. `"La_3+"`.

#Examples
```julia-repl
julia >
```
"""
const IONIC_RADIUS_CNEIGHT::Base.ImmutableDict{String,Float64} = Base.ImmutableDict(
    "La_3+" => 1.160,
    "Ce_3+" => 1.143,
    "Pr_3+" => 1.126,
    "Nd_3+" => 1.109,
    "Pm_3+" => 1.093,
    "Sm_3+" => 1.079,
    "Eu_3+" => 1.066,
    "Gd_3+" => 1.053,
    "Tb_3+" => 1.040,
    "Dy_3+" => 1.027,
    "Ho_3+" => 1.015,
    "Er_3+" => 1.004,
    "Tm_3+" => 0.994,
    "Yb_3+" => 0.985,
    "Lu_3+" => 0.977
)

#Isotopic properties (mass, radii, charge)

#C1 Chondrite (Palme and O'Neill 2016)
const CI_CHONDRITE_PO2016::Base.ImmutableDict{Symbol,Vector{Float64}} = Base.ImmutableDict(
    :H => [19700.0, 1970.0],
    :Li => [1.45, 0.145],
    :Be => [0.0219, 0.001533],
    :B => [0.775, 0.0775],
    :C => [34800.0, 3480.0],
    :N => [2950.0, 442.5],
    :O => [459000.0, 45900.0],
    :F => [58.2, 9.312],
    :Na => [4962.0, 446.58],
    :Mg => [95400.0, 3816.0],
    :Al => [8400.0, 504.0],
    :Si => [107000.0, 3210.0],
    :P => [985.0, 78.8],
    :S => [53500.0, 2675.0],
    :Cl => [698.0, 104.7],
    :K => [546.0, 49.14],
    :Ca => [9110.0, 546.6],
    :Sc => [5.81, 0.3486],
    :Ti => [447.0, 31.29],
    :V => [54.6, 3.276],
    :Cr => [2623.0, 131.15],
    :Mn => [1916.0, 114.96],
    :Fe => [186600.0, 7464],
    :Co => [513.0, 20.52],
    :Ni => [10910.0, 763.7],
    :Cu => [133.0, 18.62],
    :Zn => [309.0, 12.36],
    :Ga => [9.62, 0.5772],
    :Ge => [32.6, 2.934],
    :As => [1.74, 0.1566],
    :Se => [20.3, 1.421],
    :Br => [3.26, 0.489],
    :Rb => [2.32, 0.1856],
    :Sr => [7.79, 0.5453],
    :Y => [1.46, 0.073],
    :Zr => [3.63, 0.1815],
    :Nb => [0.283, 0.0283],
    :Mo => [0.961, 0.0961],
    :Ru => [0.69, 0.0345],
    :Rh => [0.132, 0.0066],
    :Pd => [0.56, 0.0224],
    :Ag => [0.201, 0.01809],
    :Cd => [0.674, 0.04718],
    :In => [0.0778, 0.00389],
    :Sn => [1.63, 0.2445],
    :Sb => [0.145, 0.0203],
    :Te => [2.28, 0.1596],
    :I => [0.53, 0.106],
    :Cs => [0.188, 0.01128],
    :Ba => [2.42, 0.121],
    :La => [0.2414, 0.007242],
    :Ce => [0.6194, 0.018582],
    :Pr => [0.0939, 0.002817],
    :Nd => [0.4737, 0.014211],
    :Sm => [0.1536, 0.004608],
    :Eu => [0.05883, 0.0017649],
    :Gd => [0.2069, 0.006207],
    :Tb => [0.03797, 0.0011391],
    :Dy => [0.2558, 0.007674],
    :Ho => [0.05644, 0.0016932],
    :Er => [0.1655, 0.004965],
    :Tm => [0.02609, 0.0007827],
    :Yb => [0.1687, 0.005061],
    :Lu => [0.02503, 0.0007509],
    :Hf => [0.1065, 0.003195],
    :Ta => [0.015, 0.0015],
    :W => [0.096, 0.0096],
    :Re => [0.04, 0.002],
    :Os => [0.495, 0.02475],
    :Ir => [0.469, 0.02345],
    :Pt => [0.925, 0.04625],
    :Au => [0.148, 0.01776],
    :Hg => [0.35, 0.175],
    :Tl => [0.14, 0.0154],
    :Pb => [2.62, 0.2096],
    :Bi => [0.11, 0.0099],
    :Th => [0.03, 0.0021],
    :U => [0.0081, 0.000567],
)

# Initial Sr87/86 ratios
"""
    INITIAL_Sr["key"]
Return the initial ⁸⁷Sr/⁸⁶Sr or ⁸⁶Sr/⁸⁷Sr for a given reference material and its standard error.

type 'keys(INITIAL_Sr)' for a list of available materials
"""
const INITIAL_Sr::Base.ImmutableDict{String,Vector{Float64}} = Base.ImmutableDict(
    "MDC" => [0.72607, 0.00035],
    "MDCInv" => [1.37727767295164, 0.00066391282594388300],
    "Hogsbo" => [0.72000, 0.0000515],
    "HogsboInv" => [1.38888888888888, 0.00009934413580246910],
    "MicaMg" => [0.72607, 0.00035],
    "MicaMgInv" => [1.37727767295164, 0.00066391282594388300],
    "SolarSystem" => [0.69899, 0.00005],
    "SolarSystemInv" => [1.43063563141103, 0.00010233591549314200],
    "GLO" => [0.707358, 0.00005],
    "GLOInv" => [1.41371130318735, 0.00009992898243798420],
    "LaPosta" => [0.70483, 0.00005],
    "LaPostaInv" => [1.41878183391739, 0.00010064709461270000],
    "MtDromedary" => [0.7049, 0.00005],
    "MtDromedaryInv" => [1.41864094197758, 0.00010062710611275200],
    "WilsonsProm" => [0.71, 0.00005],
    "WilsonsPromInv" => [1.40845070422535, 0.00009918666931164450],
    "NIST610" => [0.709699, 0.000009],
    "NIST610Inv" => [1.40904806122032, 0.00001786874794945860],
    "BCR2G" => [0.705003, 0.000004],
    "BCR2GInv" => [1.41843368042405, 0.00000804781642304532],
    "NIST612" => [0.7090630, 0.00001],
    "NIST612Inv" => [1.41031191868705, 0.00001988979707990760]
)
