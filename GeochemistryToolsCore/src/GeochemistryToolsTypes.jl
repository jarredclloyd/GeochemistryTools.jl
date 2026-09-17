# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd# SPDX-FileCopyrightText: Copyright © 2026Jarred C Lloyd 4o8pg2v9h@mozmail.com
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

# GeochemistryToolsCore.jl/GeochemistryToolsTypes.jl
# High level struct and type definitions for custom structs and abstract types used across the GeochemistryTools ecosystem.

export LinearRegression, OrthogonalPolynomial, Mineral, Nuclide

abstract type LinearRegression <: Any end

abstract type OrthogonalPolynomial <: LinearRegression end

abstract type Mineral <: Any end

abstract type Nuclide <: Any end

Base.getindex(LR::LinearRegression, key::Symbol) = getfield(LR, key)
Base.getindex(LR::LinearRegression, key::AbstractString) = getfield(LR, Symbol(key))
Base.setindex!(LR::LinearRegression, value, key::Symbol) = setfield!(LR, key, value)
function Base.setindex!(LR::LinearRegression, value, key::AbstractString)
    return setfield!(LR, Symbol(key), value)
end
Base.haskey(LR::LinearRegression, key::Symbol) = key ∈ fieldnames(typeof(LR))
Base.keys(LR::LinearRegression) = fieldnames(typeof(LR))
Base.values(LR::LinearRegression) = getfield.(Ref(LR), keys(LR))
Base.pairs(LR::LinearRegression) = zip(keys(LR), values(LR))


Base.getindex(mineral::Mineral, key::Symbol) = getfield(mineral, key)
Base.getindex(mineral::Mineral, key::AbstractString) = getfield(mineral, Symbol(key))
Base.setindex!(mineral::Mineral, value, key::Symbol) = setfield!(mineral, key, value)
function Base.setindex!(mineral::Mineral, value, key::AbstractString)
    return setfield!(mineral, Symbol(key), value)
end
Base.haskey(mineral::Mineral, key::Symbol) = key ∈ fieldnames(typeof(mineral))
Base.keys(mineral::Mineral) = fieldnames(typeof(mineral))
Base.values(mineral::Mineral) = getfield.(Ref(mineral), keys(mineral))
Base.pairs(mineral::Mineral) = zip(keys(mineral), values(mineral))

Base.getindex(nuclide::Nuclide, key::Symbol) = getfield(nuclide, key)
Base.getindex(nuclide::Nuclide, key::AbstractString) = getfield(nuclide, Symbol(key))
Base.setindex!(nuclide::Nuclide, value, key::Symbol) = setfield!(nuclide, key, value)
function Base.setindex!(nuclide::Nuclide, value, key::AbstractString)
    return setfield!(nuclide, Symbol(key), value)
end
Base.haskey(nuclide::Nuclide, key::Symbol) = key ∈ fieldnames(typeof(nuclide))
Base.keys(nuclide::Nuclide) = fieldnames(typeof(nuclide))
Base.values(nuclide::Nuclide) = getfield.(Ref(nuclide), keys(nuclide))
Base.pairs(nuclide::Nuclide) = zip(keys(nuclide), values(nuclide))
