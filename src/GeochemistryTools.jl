# GeochemistryTools.jl
# Functions for geochemical and geochronological data in julia
#
# Copyright © 2023 Jarred C Lloyd
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the “Software”), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
__precompile__()

module GeochemistryTools

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
@reexport using Dates
@reexport using SparseArrays
using LinearAlgebra
using HypergeometricFunctions
using SpecialFunctions
using PyCall
using Conda

import Base: getindex, setindex!
import Base.Threads: @spawn, @threads, @simd

include("formula_to_weight.jl")
include("GCTDictionaries.jl")
include("Regression/ErrorInVariablesRegression/errors_in_variable_regression.jl")
include("Geochronology/beta_minus_decay_systems.jl")
include("Geochronology/UPb.jl")
include("ICP_MS.jl")
include("raman_spectroscopy.jl")
include("profilometry.jl")
include("Minerals/mineral_formulas.jl")
include("lanthanoid_lambdas.jl")
include("Regression/ErrorInVariablesRegression/eivlr_deming.jl")
include("Regression/ErrorInVariablesRegression/eivlr_mahon.jl")
include("Regression/ErrorInVariablesRegression/eivlr_york.jl")
include("Regression/generalised_least_squares.jl")
include("Regression/measures_of_fit.jl")
include("Regression/orthogonal_polynomials.jl")
include("geometric_statistics.jl")
include("error_ellipse.jl")

function _check_equal_length(
    a::AbstractVector,
    b::AbstractVector,
    c::Union{T, Nothing}=nothing,
    d::Union{T, Nothing}=nothing) where {T<:AbstractVector}
    if c !== nothing && d !== nothing
        length(a) == length(b) == length(c) == length(d)
    elseif c !== nothing && d === nothing
        length(a) == length(b) == length(c)
    else
        length(a) == length(b)
    end
end

const pybaselines = PyNULL()

pyimport_conda("scipy", "scipy")

function __init__()
    copy!(pybaselines, pyimport_conda("pybaselines", "pybaselines"))
end

println("Hello 👋 \n If you wish to use the 'plot' functions of this package you will need to add a Makie backend (e.g. GLMakie, CairoMakie)")

end
