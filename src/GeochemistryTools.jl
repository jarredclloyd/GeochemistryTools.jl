module GeochemistryTools

using Reexport
using Base.Threads: @spawn, @threads
@reexport using CSV
@reexport using DataFrames
@reexport using Distributions
@reexport using Polynomials
@reexport using Statistics
@reexport using StatsBase
@reexport using ColorSchemes
@reexport using Glob
@reexport using GLMakie
@reexport using HypothesisTests

using PyCall
using Conda

include("FormulaToWeight.jl")
include("GCTDictionaries.jl")
include("UPb.jl")
include("EIVLinearRegression.jl")
include("RbSr.jl")
include("ICP_MS.jl")
include("RamanSpectroscopy.jl")

const pybaselines = PyNULL()

pyimport_conda("scipy", "scipy")

function __init__()
    copy!(pybaselines, pyimport_conda("pybaselines", "pybaselines"))
end

end