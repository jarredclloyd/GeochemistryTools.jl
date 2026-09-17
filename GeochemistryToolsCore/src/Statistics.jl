# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

# GeochemistryToolsCore.jl/Statistics.jl
# Statistical functions to work with uncertainty tracking and long term data

export logmean, logstdmean, weightedmean, GrandMeanResult, grandmean, geoweightedmean

function logmean(𝑥)
    @assert all(𝑥 .> 0) "values ≤ 0 are present"
    if typeof(𝑥[begin]) <: Measurement
        return mean(log.(Measurement.value.(𝑥)))
    else
        return mean(log.(𝑥))
    end
end

function logstdmean(𝑥) # experimental standard deviation of the mean (aka standard error)
    @assert all(𝑥 .> 0) "values ≤ 0 are present"
    if typeof(𝑥[begin]) <: Measurement
        return std(log.(Measurement.value.(𝑥))) / sqrt(length(𝑥))
    else
        return std(log.(𝑥)) / sqrt(length(𝑥))
    end
end

function weightedmean(𝑥)
    𝜔 = 1 ./ Measurements.uncertainty.(𝑥) .^ 2
    return sum(𝜔 .* 𝑥) / sum(𝜔)
end

function geoweightedmean(𝑥)
    @assert all(𝑥 .> 0) "values ≤ 0 are present"
    return exp(GeochemistryToolsCore.weightedmean(log.(𝑥)))
end

Base.@kwdef struct GrandMeanResult
    mean::Measurement
    tau::Measurement
    chi2::Float64
    dof::Int
    outlier_inds::Union{Missing,Vector{Int}}
    reduced_chi2_in::Float64
    n_total::Int
    n_fit::Int
end

Base.getindex(gm::GrandMeanResult, key::Symbol) = getfield(gm, key)
function Base.getindex(gm::GrandMeanResult, key::AbstractString)
    return getfield(gm, Symbol(key))
end
function Base.haskey(gm::GrandMeanResult, key::Symbol)
    return key ∈ fieldnames(typeof(gm))
end
Base.keys(gm::GrandMeanResult) = fieldnames(typeof(gm))
Base.values(gm::GrandMeanResult) = getfield.(Ref(gm), keys(gm))
Base.pairs(gm::GrandMeanResult) = zip(keys(gm), values(gm))

function Base.show(io::IOContext, gm::GrandMeanResult)
    if get(io, :compact, false)::Bool
        return print(io, gm.mean, " (τ = $(gm.tau * gm.mean.val))")
    else
        pretty_table(
            hcat(
                gm.mean,
                gm.tau * 100,
                gm.n_total,
                gm.n_fit,
                string(gm.outlier_inds),
                gm.reduced_chi2_in,
                gm.chi2 ./ gm.dof
            );
            title = "Grand mean summary",
            column_labels = [
                "Geometric Mean",
                "Excess Variance (𝜏) 1s%",
                "𝑛 total",
                "𝑛 fit",
                "outlier indices",
                "χ²ᵣ (input)",
                "χ²ᵣ (including 𝜏)"
            ]
        )
    end
end

"""
    grandmean(x::Vector{Measurement{Float64}}; kwargs...)

    Compute the grand mean of independent values using a variance components model that estimates an excess variance term 𝜏² such that the reduced χ² (MSWD) approaches 1.

    # Keywords
    - `rm_outlier::Bool=true`: Enable automated outlier detection using modified z-score
    - `outlier_tol::Real=2`: Z-score tolerance (in standard deviations)
"""
function grandmean(
    x::Vector{Measurement{Float64}};
    rm_outlier::Bool = true,
    outlier_tol::Real = 2
)
    𝑛ₜ = length(x)
    outliers::BitVector = falses(𝑛ₜ)
    if rm_outlier
        MAD = mad(x; normalize = false)
        if !iszero(MAD)
            𝐳 = 0.6745 .* (Measurements.value.(x) .- median(x)) ./ MAD
            outliers .= abs.(𝐳) .> outlier_tol
        end
    end
    out_mask = .!outliers
    vals = Measurements.value.(x[out_mask])
    uncs = Measurements.uncertainty.(x[out_mask])
    𝜈 = length(vals) - 1
    𝜇̄::Measurement = GeochemistryToolsCore.weightedmean(measurement.(vals, uncs))
    𝛚 = 1 ./ (uncs .^ 2)
    χ² = sum(𝛚 .* (vals .- 𝜇̄.val) .^ 2)
    χ²ᵣ::Float64 = χ² / 𝜈
    𝐱 = log.(vals)
    𝐮 = uncs ./ vals
    𝛚 = 1 ./ (𝐮 .^ 2)
    𝜇̄ = sum(𝛚 .* 𝐱) / sum(𝛚)
    χ² = sum(𝛚 .* (𝐱 .- 𝜇̄.val) .^ 2)
    𝜏² = max(0, (χ² - 𝜈) / (sum(𝛚) - sum(𝛚 .^ 2) / sum(𝛚)))
    if χ² ≤ 𝜈
        𝜏² = 0.0
        𝜏²ₛ = 0.0
    else
        𝜏², misfit_deriv = _refine_𝜏²(𝐱, 𝐮, χ² / 𝜈; 𝜏² = 𝜏²)
        𝜏²ₛ = sqrt(2𝜈) / abs(misfit_deriv)
    end

    𝛚ₐ = 1.0 ./ (𝐮 .^ 2 .+ 𝜏²)
    μ̂ₐ = sum(𝛚ₐ .* 𝐱) / sum(𝛚ₐ)
    χ² = sum(𝛚ₐ .* (𝐱 .- μ̂ₐ) .^ 2)

    𝐮ₐ = sqrt.(𝐮 .^ 2 .+ 𝜏²)
    𝐱𝐬ₐ = measurement.(𝐱, 𝐮ₐ)
    𝜇̄ = exp(GeochemistryToolsCore.weightedmean(𝐱𝐬ₐ))
    𝜏 = sqrt(𝜏²)
    𝜏ₛ = iszero(𝜏) ? 0 : 𝜏²ₛ / (2 * 𝜏)

    return GrandMeanResult(;
        mean = 𝜇̄,
        tau = 𝜏 ± 𝜏ₛ,
        chi2 = χ²,
        dof = 𝜈,
        outlier_inds = findall(outliers),
        reduced_chi2_in = χ²ᵣ,
        n_total = 𝑛ₜ,
        n_fit = 𝜈 + 1
    )
end

function _refine_𝜏²(x, u, χ²ᵣ::Float64 = Inf; 𝜏²::Float64 = 0.0, maxiter::Int = 100)
    𝜈::Int = length(x) - 1
    iter::Int = 0
    misfit_deriv::Float64 = 0
    tol::Float64 = sqrt(eps(Float64))

    while abs(χ²ᵣ - 1.0) > tol
        iter += 1
        𝐮ₐ = u .^ 2 .+ 𝜏²
        𝛚 = 1.0 ./ 𝐮ₐ

        𝜇̄ = sum(𝛚 .* x) / sum(𝛚)
        𝐫 = x .- 𝜇̄

        χ²::Float64 = sum(𝛚 .* 𝐫 .^ 2)
        χ²ᵣ = χ² / 𝜈
        χ²_misfit::Float64 = χ² - 𝜈

        # derivative
        Σ𝛚 = sum(𝛚)
        weighted_lev_grad = sum(𝐫 ./ 𝐮ₐ .^ 2)
        δ𝜇̄ = -weighted_lev_grad / Σ𝛚

        misfit_deriv = -sum(𝐫 .^ 2 ./ 𝐮ₐ .^ 2) + 2 * sum(𝐫 ./ 𝐮ₐ .* δ𝜇̄)
        if abs(misfit_deriv) < eps()
            break
        end
        𝜏²_new = 𝜏² - χ²_misfit / misfit_deriv

        if !isfinite(𝜏²_new) || 𝜏²_new < 0
            𝜏²_new = 𝜏² * 0.5
        end

        𝜏² = 𝜏²_new

        if iter ≥ maxiter
            break
        end

    end
    return max(𝜏², 0.0), misfit_deriv
end
