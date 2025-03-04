# SPDX-FileCopyrightText: 2025 Jarred Lloyd (https://github.com/jarredclloyd)
# SPDX-License-Identifier: MIT
# Created: Tuesday, March 4, 2025 @ 09:24:46 AM
# Edited: Tuesday, March 4, 2025 @ 09:24:50 AM

#= Preamble GeochemistryTools.jl - ICP-MS corrections module
This module provides functions for correcting ICP-MS data.
=#


export correct_interference
"""
    correct_interference(args...; kwargs...)

    Perform interference correction on ICP-MS data with user specified parameters.

"""

function correct_interference(
    data::DataFrame,
    distorted_signal::Union{Symbol, Integer, AbstractString, AbstractChar}, interference::Union{Symbol, Integer, AbstractString, AbstractChar},
    isotope_ratio::AbstractFloat,
    mass_bias::AbstractFloat)

    transform!(
        data,
        Cols(distorted_signal, interference) =>
            ByRow(
                (ds, int) ->
                ds - int * isotope_ratio * mass_bias) => distorted_signal * "_corr"
    )
end
