#=
This jl file contains dictionaries and other reference information for that are used in "GeochemistryTools.jl"
=#

export AVOGADRO,
    INITIAL_Sr,
    CI_CHONDRITE_PO2014,
    PAAS_P2012,
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
    "La_3+" => 116.0,
    "Ce_3+" => 114.3,
    "Pr_3+" => 112.6,
    "Nd_3+" => 110.9,
    "Pm_3+" => 109.3,
    "Sm_3+" => 107.9,
    "Eu_3+" => 106.6,
    "Gd_3+" => 105.3,
    "Tb_3+" => 104.0,
    "Dy_3+" => 102.7,
    "Ho_3+" => 101.5,
    "Er_3+" => 100.4,
    "Tm_3+" => 099.4,
    "Yb_3+" => 098.5,
    "Lu_3+" => 097.7
)

#Isotopic properties (mass, radii, charge)

# C1 Chondrite (Palme and O'Neill 2014)
const CI_CHONDRITE_PO2014::Base.ImmutableDict{Symbol,Vector{Float64}} = Base.ImmutableDict(
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

# Primitive Mantle/Bulk Silicate Earth (Palme and O'Neill 2014)
const PM_PO2014::Base.ImmutableDict{Symbol,Vector{Float64}} = Base.ImmutableDict(
    :H => [120.0, 24.0],
    :Li => [1.6, 0.32],
    :Be => [0.062, 0.0062],
    :B => [0.26, 0.104],
    :C => [100.0, 50.0],
    :N => [2.0, 1.0],
    :O => [443300.0, 8866.0],
    :F => [25.0, 10.0],
    :Na => [2590.0, 129.5],
    :Mg => [221700.0, 2217.0],
    :Al => [23800.0, 1904.0],
    :Si => [212200.0, 2122.0],
    :P => [87.0, 13.05],
    :S => [200.0, 80.0],
    :Cl => [30.0, 12.0],
    :K => [260.0, 39.0],
    :Ca => [26100.0, 2088.0],
    :Sc => [16.4, 1.64],
    :Ti => [1265.0, 126.5],
    :V => [86.0, 4.3],
    :Cr => [2520.0, 252.0],
    :Mn => [1050.0, 105.0],
    :Fe => [63000.0, 630.0],
    :Co => [102.0, 5.1],
    :Ni => [1860.0, 93.0],
    :Cu => [20.0, 10.0],
    :Zn => [53.5, 2.675],
    :Ga => [4.4, 0.22],
    :Ge => [1.2, 0.24],
    :As => [0.068, 0.0204],
    :Se => [0.076, 0.038],
    :Br => [0.075, 0.0375],
    :Rb => [0.605, 0.0605],
    :Sr => [22.0, 1.1],
    :Y => [4.13, 0.413],
    :Zr => [10.3, 1.03],
    :Nb => [0.595, 0.119],
    :Mo => [0.047, 0.0188],
    :Ru => [0.0074, 0.00148],
    :Rh => [0.0012, 0.00024],
    :Pd => [0.0071, 0.00142],
    :Ag => [0.006, 0.003],
    :Cd => [0.035, 0.007],
    :In => [0.018, 0.0036],
    :Sn => [0.14, 0.042],
    :Sb => [0.0054, 0.00216],
    :Te => [0.009, 0.0045],
    :I => [0.007, 0.0035],
    :Cs => [0.018, 0.009],
    :Ba => [6.85, 1.0275],
    :La => [0.6832, 0.06832],
    :Ce => [1.7529, 0.17529],
    :Pr => [0.2657, 0.039855],
    :Nd => [1.341, 0.1341],
    :Sm => [0.4347, 0.04347],
    :Eu => [0.1665, 0.01665],
    :Gd => [0.5855, 0.029275],
    :Tb => [0.1075, 0.016125],
    :Dy => [0.7239, 0.07239],
    :Ho => [0.1597, 0.023955],
    :Er => [0.4684, 0.04684],
    :Tm => [0.07383, 0.0110745],
    :Yb => [0.4774, 0.04774],
    :Lu => [0.07083, 0.0106245],
    :Hf => [0.3014, 0.03014],
    :Ta => [0.043, 0.00215],
    :W => [0.012, 0.0036],
    :Re => [0.00035, 0.00007],
    :Os => [0.0039, 0.000585],
    :Ir => [0.0035, 0.00035],
    :Pt => [0.0076, 0.00152],
    :Au => [0.0017, 0.00051],
    :Hg => [0.006, 0.003],
    :Tl => [0.0041, 0.001025],
    :Pb => [0.185, 0.0185],
    :Bi => [0.003, 0.0015],
    :Th => [0.0849, 0.012735],
    :U => [0.0229, 0.003435],
)

# Post-Archaean Australian Shale (Pourmand et al. 2012)
const PAAS_P2012::Base.ImmutableDict{Symbol,Vector{Float64}} = Base.ImmutableDict(
    :Sc	=> [15.89, 2.2587],
    :Y	=> [27.31, 5.3647],
    :La	=> [44.56, 7.6188],
    :Ce	=> [88.25, 15.7763],
    :Pr	=> [10.15, 1.8385],
    :Nd	=> [37.32, 7.0715],
    :Sm	=> [6.884, 1.4121],
    :Eu	=> [1.215, 0.2819],
    :Gd	=> [6.043, 1.304],
    :Tb	=> [0.8914, 0.1884],
    :Dy	=> [5.325, 1.1156],
    :Ho	=> [1.053, 0.2165],
    :Er	=> [3.075, 0.6184],
    :Tm	=> [0.4510, 0.0894],
    :Yb	=> [3.012, 0.5821],
    :Lu	=> [0.4386, 0.0821],
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
