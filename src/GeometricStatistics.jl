# SPDX-FileCopyrightText: 2024 Jarred Lloyd (https://github.com/jarredclloyd)
# SPDX-License-Identifier: MIT

#= Preamble
Last updated: 2024-10-17

This source file contains functions to compute geometric means and variances based on
Habib (2012).

Habib, EAE (2012) 'Geometric Mean for Negative and Zero Values',
International Journal of Research and Reviews in Applied Sciences, 11(3)
=#

# function exports
export geomean_zeros,
    geovar_zeros, geostd_zeros, geosem_zeros, deltalognormal, pseudolog, pseudolog_undo

"""
    geomean_zeros(x::AbstractVector)

    Computes the geometric mean for x ≥ 0.
"""
function geomean_zeros(x::AbstractVector)
    N = length(x)
    N2 = count(x .> 0)
    if N2 > 0
        G₊ = geomean(x[x[:] .> 0, :])
        G = (N2 / N) * G₊
    else
        G = 0
    end
    return G

end

"""
    geovar_zeros(x::AbstractVector)

    Computes the variance of the geometric mean for x ≥ 0
"""
function geovar_zeros(x::AbstractVector)
    N = length(x)
    N2 = count(x .> 0)
    if N2 > 1
        var₊ = exp(var(log.(x[x[:] .> 0, :])))
        varG = (N2 / N) * var₊
    else
        varG = 0
    end
    return varG
end

function geostd_zeros(x::AbstractVector)
    return std(log.(x[x[:] .> 0, :]))
end

function geosem_zeros(x::AbstractVector)
    N = length(x)
    N2 = count(x .> 0)
    if N2 > 1
        sem₊ = exp(sem(log.(x[x[:] .> 0, :])))
        semG = (N2 / N) * sem₊
    else
        semG = 0
    end
    return semG
end

"""
    tri_geomean(x::AbstractVector)

    Computes the geometric mean for x ∈ R. Only defined for odd values of 𝑁ₜ(negative)
    and 𝑁ₜ(total)
"""
function tri_geomean(x::AbstractVector)
    𝑁ₜ = length(x)
    𝑁₁ = count(x .< 0)
    𝑁₂ = count(x .> 0)
    # if 𝑁₁ > 0 && isodd(𝑁₁) !== true && isodd(𝑁ₜ) !== true
    #     throw(
    #         DomainError(
    #             "Trigeometric mean is not defined for odd total N and odd N1 (x < 0)",
    #         ),
    #     )
    # end
    𝑁₁ > 0 ? G₋ = -(geomean(abs.(x[x[:] .< 0, :]))) : G₋ = 0
    𝑁₂ > 0 ? G₊ = geomean(x[x[:] .> 0, :]) : G₊ = 0
    return (𝑁₁ * G₋ + 𝑁₂ * G₊) / 𝑁ₜ
end

function tri_geovar(x::AbstractVector)
    𝑁ₜ = length(x)
    𝑁₁ = count(x .< 0)
    𝑁₂ = count(x .> 0)
    if 𝑁₁ > 0 && isodd(𝑁₁) !== true && isodd(𝑁ₜ) !== true
        throw(
            DomainError(
                "Trigeometric mean is not defined for odd total N and odd N1 (x < 0)",
            ),
        )
    end
    if 𝑁₁ > 0
        𝑁₁𝑁 = 𝑁₁ / 𝑁ₜ
        G₋ = geomean(abs.(x[x[:] .< 0, :]))
        varlogG₋ = -(var(log.(abs.(x[x[:] .< 0, :]))))
        𝐸G₋ = _expectation_g_negative(G₋, varlogG₋, 𝑁₁, 𝑁ₜ)
        varG₋ = (𝑁₁𝑁)^2 * (varlogG₋ / 𝑁₁) * (G₋)^2
    else
        𝑁₁𝑁 = 0.0
        G₋ = 0.0
        varlogG₋ = 0.0
        𝐸G₋ = 0
        varG₋ = 0.0
    end
    if 𝑁₂ > 0
        𝑁₂𝑁 = 𝑁₂ / 𝑁ₜ
        G₊ = geomean(x[x[:] .> 0, :])
        varlogG₊ = var(log.(abs.(x[x[:] .> 0, :])))
        𝐸G₊ = _expectation_g_positive(G₊, varlogG₊, 𝑁₂, 𝑁ₜ)
        varG₊ = (𝑁₂𝑁)^2 * (varlogG₊ / 𝑁₂) * (G₊)^2
    else
        𝑁₂𝑁 = 0.0
        G₊ = 0.0
        varlogG₊ = 0.0
        𝐸G₊ = 0.0
        varG₊ = 0.0
    end
    G = (𝑁₁ * G₋ + 𝑁₂ * G₊) / 𝑁ₜ
    𝐸G = _expectation_g_weighted(G, varlogG₋, varlogG₊, 𝑁₁, 𝑁₂, 𝑁₁𝑁, 𝑁₂𝑁)
    if 𝑁₁ == 0 || 𝑁₂ == 0
        covG = 0
    else
        covG = 𝐸G - 𝐸G₋ * 𝐸G₊
    end
    return (𝑁₁𝑁)^2 * varG₋ + (𝑁₂𝑁)^2 * varG₊ - ((2 * 𝑁₁ * 𝑁₂) / 𝑁ₜ^2) * covG
end

function _expectation_g_positive(
    G₊::AbstractFloat,
    varlogG₊::AbstractFloat,
    𝑁₂::Integer,
    𝑁ₜ::Integer,
)
    𝑁₂𝑁 = 𝑁₂ / 𝑁ₜ
    return G₊ + (𝑁₂𝑁)^2 * (varlogG₊ / (2 * 𝑁₂)) * G₊
end

function _expectation_g_negative(
    G₋::AbstractFloat,
    varlogG₋::AbstractFloat,
    𝑁₁::Integer,
    𝑁ₜ::Integer,
)
    𝑁₁𝑁 = 𝑁₁ / 𝑁ₜ
    return G₋ + (𝑁₁𝑁)^2 * (varlogG₋ / (2 * 𝑁₁)) * G₋
end

function _expectation_g_weighted(
    G::AbstractFloat,
    varlogG₋::AbstractFloat,
    varlogG₊::AbstractFloat,
    𝑁₁::Integer,
    𝑁₂::Integer,
    𝑁₁𝑁::AbstractFloat,
    𝑁₂𝑁::AbstractFloat,
)
    return G + (G / 2) * ((𝑁₁𝑁)^2 * (varlogG₋ / 𝑁₁) + (𝑁₂𝑁)^2 * (varlogG₊ / 𝑁₂))
end

function _geomean_zeros_cruz(x::AbstractVector, ϵ::AbstractFloat = 1e-5)
    G₊ = geomean(x[x[:] .> 0, :])
    δmin = 0
    δmax = G₊ - minimum(x[x[:] .> 0, :])
    δ = (δmin + δmax) / 2
    ϵ = ϵ * G₊
    auxExp = geomean(x[x[:] .> 0, :] .+ δ) - δ
    while (auxExp - G₊) > ϵ
        auxExp < G₊ ? δmin = δ : δmax = δ
        δ = (δmin + δmax) / 2
        auxExp = geomean(x[x[:] .> 0, :] .+ δ) - δ
    end
    G = geomean(x .+ δ) - δ
    return (G, δ)
end

"""
    deltalognormal(x::AbstractVector)

    Compute the mean and variance of a delta-lognormal distribution from 'x'
"""
function deltalognormal(x::AbstractVector{T}) where {T<:Real}
    𝑁 = length(x[isfinite.(x)])
    𝑁z = length(x[isfinite.(x) .&& x .== 0])
    if 𝑁z == 𝑁
        γ = 0
        δ = 0
        return γ, δ
    else
        𝜃 = 𝑁z / 𝑁
        if 𝑁 - 𝑁z == 1
            μ, σ = log(x[isfinite.(x) .&& x .> 0][1]), log(1.0)
            γ = (1 - 𝜃) * exp(μ + σ^2 / 2)
            δ = (1 - 𝜃) * exp(2 * μ + σ^2) * exp(σ^2 - (1 - 𝜃))
        else
            dlog = Distributions.fit(LogNormal, x[isfinite.(x) .&& x .> 0])
            γ = (1 - 𝜃) * exp(dlog.μ + dlog.σ^2 / 2)
            δ = (1 - 𝜃) * exp(2 * dlog.μ + dlog.σ^2) * exp(dlog.σ^2 - (1 - 𝜃))
        end
        return γ, δ
    end
end

function pseudolog(x, σ = 0.05, base = ℯ)
    return asinh(x / (2 * σ) / log(base))
end

function pseudolog_undo(x, σ = 0.05, base = ℯ)
    return 2 * σ * sinh(x * log(base))
end
