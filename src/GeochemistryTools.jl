# SPDX-FileCopyrightText: 2024 Jarred Lloyd (https://github.com/jarredclloyd)
# SPDX-License-Identifier: MIT

# GeochemistryTools.jl
# Functions for geochemical and geochronological data in julia

__precompile__()

module GeochemistryTools

using Reexport
@reexport using CSV
@reexport using DataFrames
using Distributions
using Polynomials
using Statistics
using StatsBase
@reexport using ColorSchemes
using Glob
using HypothesisTests
@reexport using Dates
using SparseArrays
using LinearAlgebra
using GenericLinearAlgebra
using MultiFloats
using SpecialFunctions
using PyCall
using Conda
using PolynomialRoots
@reexport using PeriodicTable
@reexport using IsotopeTable
@reexport using Unitful

import Base: getindex, setindex!
import Base.Threads: @spawn, @threads, @simd

include("GCTDictionaries.jl")
include("Regression/ErrorInVariablesRegression/ErrorsInVariablesRegression.jl")
include("Geochronology/beta_minus_decay_systems.jl")
include("Geochronology/UPb.jl")
include.(filter(contains(r".jl$"), readdir(joinpath(normpath(@__DIR__),"ICPMS/"); join=true)))
include("RamanSpectroscopy.jl")
include("Profilometry.jl")
include("Minerals/minerals.jl")
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
include("BarycentricConversions.jl")
include.(filter(contains(r".jl$"), readdir(joinpath(normpath(@__DIR__),"ChemistryTools/"); join=true)))

# To enable use of MultiFloats in Julia 1.11 until MultiFloats v3.0 is released
@inline Base.precision(::Type{MultiFloat{T,N}}) where {T,N} =
           N * precision(T) + (N - 1) # implicit bits of precision between limbs
@inline    Base.round(x::MultiFloat{T,N} where {T,N}, r::RoundingMode) = Base.round(big.(x), r)
@inline    Base.trunc(x::MultiFloat{T,N} where {T,N}, r::RoundingMode{:ToZero}) = Base.trunc(big.(x), r)

function _check_equal_length(
    a::AbstractVector,
    b::AbstractVector,
    c::Union{AbstractVector,Nothing} = nothing,
    d::Union{AbstractVector,Nothing} = nothing,
)
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


# switch to PythonCall and CondaPkg when open_ssl is updated, or write baselines into native julia
#=
const pybaselines = Ref{Py}()
function __init__
    pybaselines = pyimport("pybaselines)
end
=#
const pybaselines = PyNULL()

function __init__()
    copy!(pybaselines, pyimport_conda("pybaselines", "pybaselines"))
    MultiFloats.use_bigfloat_transcendentals()
    println(
        "👋 Thanks for using GeochemistryTools \n If you wish to use the 'plot' functions of this package you will need to add a Makie backend (e.g. GLMakie, CairoMakie)",
    )
end

end
