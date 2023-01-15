module GeochemistryTools

using Base.Threads: @spawn, @threads
using Reexport
@reexport using CSV
@reexport using DataFrames
@reexport using Distributions
@reexport using Polynomials
@reexport using Statistics
@reexport using StatsBase
@reexport using ColorSchemes
@reexport using Glob
@reexport using HypothesisTests
@reexport using GLMakie
@reexport using Dates

using PyCall
using Conda

include("FormulaToWeight.jl")
include("GCTDictionaries.jl")
include("UPb.jl")
include("EIVLinearRegression.jl")
include("RbSr.jl")
include("ICP_MS.jl")
include("RamanSpectroscopy.jl")
include("Profilometer.jl")

const pybaselines = PyNULL()

pyimport_conda("scipy", "scipy")

function __init__()
    copy!(pybaselines, pyimport_conda("pybaselines", "pybaselines"))
end

println("Hello 👋 \n If you wish to use the 'plot' functions of this package you will need to add a Makie backend 
(e.g. GLMakie, CairoMakie)")

end