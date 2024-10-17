# SPDX-FileCopyrightText: 2024 Jarred Lloyd (https://github.com/jarredclloyd)
# SPDX-License-Identifier: MIT

# GeochemistryTools.jl
# Functions for geochemical and geochronological data in julia

__precompile__()

module GeochemistryTools

using Reexport
using CSV
using DataFrames
using Distributions
using Polynomials
using Statistics
using StatsBase
using ColorSchemes
using Glob
using HypothesisTests
using Dates
using SparseArrays
using LinearAlgebra
using MultiFloats
using SpecialFunctions
using PyCall
using Conda
using PolynomialRoots

import Base: getindex, setindex!
import Base.Threads: @spawn, @threads, @simd

include("FormulaWeight.jl")
include("GCTDictionaries.jl")
include("Regression/ErrorInVariablesRegression/ErrorsInVariablesRegression.jl")
include("Geochronology/beta_minus_decay_systems.jl")
include("Geochronology/UPb.jl")
include("ICP_MS.jl")
include("RamanSpectroscopy.jl")
include("Profilometry.jl")
include("Minerals/mineral_formulas.jl")
include("LanthanoidLambdas.jl")
include("Regression/ErrorInVariablesRegression/EIVLRDeming.jl")
include("Regression/ErrorInVariablesRegression/EIVLRMahon.jl")
include("Regression/ErrorInVariablesRegression/EIVLRYork.jl")
include("Regression/GeneralisedLeastSquares.jl")
include("Regression/MeasuresOfFit.jl")
include("Regression/OrthogonalPolynomials.jl")
include("GeometricStatistics.jl")
include("ErrorEllipse.jl")
include("DateTimeParser.jl")

function _check_equal_length(
    a::AbstractVector,
    b::AbstractVector,
    c::Union{T,Nothing} = nothing,
    d::Union{T,Nothing} = nothing,
) where {T<:AbstractVector}
    if c !== nothing && d !== nothing
        length(a) == length(b) == length(c) == length(d)
    elseif c !== nothing && d === nothing
        length(a) == length(b) == length(c)
    else
        length(a) == length(b)
    end
end

function timeout(f, args, seconds, fail)
    tsk = @task f(args...)
    schedule(tsk)
    Timer(seconds) do timer
        return istaskdone(tsk) || Base.throwto(tsk, InterruptException())
    end
    try
        fetch(tsk)
    catch _
        fail
    end
end

const pybaselines = PyNULL()

pyimport_conda("scipy", "scipy")

function __init__()
    copy!(pybaselines, pyimport_conda("pybaselines", "pybaselines"))
    MultiFloats.use_bigfloat_transcendentals()
    println(
        "👋 Thanks for using GeochemistryTools \n If you wish to use the 'plot' functions of this package you will need to add a Makie backend (e.g. GLMakie, CairoMakie)",
    )
end

end
