#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-10-05
Edited: 2023-10-05

This source file contains functions to compute geometric means and variances based on
Habib (2012).

Habib, EAE (2012) 'Geometric Mean for Negative and Zero Values',
International Journal of Research and Reviews in Applied Sciences, 11(3),
www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf

=#
# function exports
export geomean_zeros, geovar_zeros

"""
    geomean_zeros(a, ϵ=1e-5; algorithm ="habib")

    Computes the geometric mean for x ≥ 0.
        Two algorithms are available: `"habib"` or `"cruz_kreft"`
"""
function geomean_zeros(a::AbstractVector)
    N = length(a)
    N2 = count(a .> 0)
    if N2 > 0
        G₊ = geomean(a[a[:] .> 0, :])
        G = (N2 / N) * G₊
    else
        G = 0
    end
    return G

end

"""
    geovar_zeros(a::AbstractVector)

    Computes the variance of the geometric mean for x ≥ 0
"""
function geovar_zeros(a::AbstractVector)
    N = length(a)
    N2 = count(a .> 0)
    if N2 > 0
        G₊ = geomean(a[a[:] .> 0, :])
        var₊ = var(log.(a[a[:] .> 0, :]))
        varG = (N2 / N)^2 * (var₊ / N2) * G₊^2
    else
        varG = 0
    end
    return varG
end

"""
    tri_geomean(a::AbstractVector)

    Computes the geometric mean for x ∈ R. Only defined for odd values of 𝑁ₜ(negative)
    and 𝑁ₜ(total)
"""
function tri_geomean(a::AbstractVector)
    𝑁ₜ = length(a)
    𝑁₁ = count(a .< 0)
    𝑁₂ = count(a .> 0)
    # if 𝑁₁ > 0 && isodd(𝑁₁) !== true && isodd(𝑁ₜ) !== true
    #     throw(
    #         DomainError(
    #             "Trigeometric mean is not defined for odd total N and odd N1 (x < 0)",
    #         ),
    #     )
    # end
    𝑁₁ > 0 ? G₋ = -(geomean(abs.(a[a[:] .< 0, :]))) : G₋ = 0
    𝑁₂ > 0 ? G₊ = geomean(a[a[:] .> 0, :]) : G₊ = 0
    return (𝑁₁ * G₋ + 𝑁₂ * G₊) / 𝑁ₜ
end

function tri_geovar(a::AbstractVector)
    𝑁ₜ = length(a)
    𝑁₁ = count(a .< 0)
    𝑁₂ = count(a .> 0)
    if 𝑁₁ > 0 && isodd(𝑁₁) !== true && isodd(𝑁ₜ) !== true
        throw(
            DomainError(
                "Trigeometric mean is not defined for odd total N and odd N1 (x < 0)",
            ),
        )
    end
    if 𝑁₁ > 0
        𝑁₁𝑁 = 𝑁₁ / 𝑁ₜ
        G₋ = geomean(abs.(a[a[:] .< 0, :]))
        varlogG₋ = -(var(log.(abs.(a[a[:] .< 0, :]))))
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
        G₊ = geomean(a[a[:] .> 0, :])
        varlogG₊ = var(log.(abs.(a[a[:] .> 0, :])))
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

function _geomean_zeros_cruz(a::AbstractVector, ϵ::AbstractFloat = 1e-5)
    G₊ = geomean(a[a[:] .> 0, :])
    δmin = 0
    δmax = G₊ - minimum(a[a[:] .> 0, :])
    δ = (δmin + δmax) / 2
    ϵ = ϵ * G₊
    auxExp = geomean(a[a[:] .> 0, :] .+ δ) - δ
    while (auxExp - G₊) > ϵ
        auxExp < G₊ ? δmin = δ : δmax = δ
        δ = (δmin + δmax) / 2
        auxExp = geomean(a[a[:] .> 0, :] .+ δ) - δ
    end
    G = geomean(a .+ δ) - δ
    return (G, δ)
end
