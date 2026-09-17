# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

# GeochemistryToolsCore.jl
# Common components for the GeochemistryTools ecosystem in Julia

module GeochemistryToolsCore

using Base.Threads: @spawn, @threads, @simd
using Downloads
using TOML
using YAML
using Reexport
using IsotopeTable
using Mendeleev
@reexport begin
    using IsotopeTable: isotopes, Isotope
    using Mendeleev: elements, chem_elements, Element
    using CSV
    using Dates
    using DataFrames
    using Distributions
    using Glob
    using HypothesisTests
    using LinearAlgebra
    using Measurements
    using MultiFloats
    using PrettyTables
    using SparseArrays
    using StaticArrays
    using Statistics
    using StatsBase
    using UUIDs
    using Unitful
end

const _file_lookup_dir = normpath(homedir(), ".julia/GeochemistryTools")
const _path_radionuclides::String = normpath(_file_lookup_dir, "Radionuclides.yaml")
const _path_refvals::String = normpath(_file_lookup_dir, "GeochemicalReferenceValues.yaml")
const _path_refmats::String =
    normpath(_file_lookup_dir, "GeochemistryReferenceMaterials.yaml")

include("GeochemistryToolsTypes.jl")
include("UnicodeMap.jl")
include("DateTimeParser.jl")
include("ReferenceValueConstants.jl")
include("Transformations.jl")
include("GeometricStatistics.jl")
include("AccessorFunctions.jl")
include("SmoothingFunctions.jl")
include("MeasuresOfFit.jl")
include("UncertaintyContours.jl")
include("UncertaintyTracking.jl")
include("Statistics.jl")

Unitful._basefactors(@__MODULE__)
"""
GeochemistryToolsCore.aʸʳ
\nThe annum, a unit of time, defined as 31556925.445s by IUPAC-IUGS
\nDimension: [`Unitful.𝐓`](@ref).
\nSee Also: [`Unitful.s`](@ref).
"""
@unit aʸʳ "a" annum 3.1556925445e7 * Unitful.s true true
"""
GeochemistryToolsCore.y

\nThe tropical year, a unit of time, defined as 31556926
\nDimension: [`Unitful.𝐓`](@ref).
\nSee Also: [`Unitful.s`](@ref).
"""
@unit y "y" year_tropical 31556926 * Unitful.s true true

function __init__()
    Unitful.register(@__MODULE__)

    if ccall(:jl_generating_output, Cint, ()) == 0
        if !ispath(_file_lookup_dir)
            mkpath(_file_lookup_dir)
        end

        _check_ref_files()

        println("Constructing `RADIONUCLIDES` with file in $(GeochemistryToolsCore._file_lookup_dir)",)
        RADIONUCLIDES[] = _construct_radionuclides(_path_radionuclides)

        println("Constructing `REFERENCE_MATERIALS` with file in $(GeochemistryToolsCore._file_lookup_dir)",)
        REFERENCE_MATERIALS[] = _construct_reference_materials(_path_refmats)

        println("Constructing `GEOCHEM_REFVALS` with file in $(GeochemistryToolsCore._file_lookup_dir)",)
        GEOCHEM_REFVALS[] = _construct_reference_values(_path_refvals)

        return println("Thanks for using GeochemistryTools.jl")
    end
end

export delete_geochemistry_tools_folder!
"""
    delete_geochemistry_tools_folder!()

    Function to delete $_file_lookup_dir and all files within it.
    Use when removing GeochemistryTools or when you want to reset the reference files and folder.
    Will prompt user for confirmation.

    See also: `update_reference_files!()`
"""
function delete_geochemistry_tools_folder!()
    println("This will remove the folder: $_file_lookup_dir and all files within it.\nDo you want to proceed [Y/N]",)
    answer = readline()
    if occursin("y", lowercase(answer))
        rm(_file_lookup_dir; recursive = true)
    else
        println("removal aborted")
    end
end

end
