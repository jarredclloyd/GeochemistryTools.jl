module GeochemistryTools

using Statistics, StatsBase, CSV, DataFrames, CairoMakie, GLMakie
using Base.Threads: @spawn, @threads

include("FormulaToWeight.jl")
include("GCTDictionaries.jl")
include("UPb.jl")

end
