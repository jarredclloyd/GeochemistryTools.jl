module GeochemistryTools

using Statistics, StatsBase, CSV, DataFrames, GLM, Distributions
using Base.Threads: @spawn, @threads

include("FormulaToWeight.jl")
include("GCTDictionaries.jl")
include("UPb.jl")
include("EIVLinearRegression.jl")
include("RbSr.jl")

end
