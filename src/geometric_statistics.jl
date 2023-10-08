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
    geomean_zeros(a::AbstractVector)

    Computes the geometric mean for x ≥ 0
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
    var₊ = var(log.(a[a[:] .> 0, :]))
    return (N2 / N)^2 * var₊
end

"""
    tri_geomean(a::AbstractVector)

    Computes the geometric mean for x ∈ R. Only defined for odd values of 𝑁(negative)
    and 𝑁(total)
"""
function tri_geomean(a::AbstractVector)
    𝑁 = length(a)
    𝑁₁ = count(a .< 0)
    𝑁₂ = count(a .> 0)
    if 𝑁₁ > 0 && isodd(𝑁₁) !== true && isodd(𝑁) !== true
        throw(
            DomainError(
                "Trigeometric mean is not defined for odd total N and odd N1 (x < 0)",
            ),
        )
    end
    if 𝑁₁ > 0
        G₋ = -(geomean(abs.(a[a[:] .< 0, :])))
    else
        G₋ = 0
    end
    if 𝑁₂ > 0
        G₊ = geomean(a[a[:] .> 0, :])
    else
        G₊ = 0
    end
    return (𝑁₁ * G₋ + 𝑁₂ * G₊) / 𝑁
end

function tri_geovar(a::AbstractVector)
    𝑁 = length(a)
    𝑁₁ = count(a .< 0)
    𝑁₂ = count(a .> 0)
    if 𝑁₁ > 0 && isodd(𝑁₁) !== true && isodd(𝑁) !== true
        throw(
            DomainError(
                "Trigeometric mean is not defined for odd total N and odd N1 (x < 0)",
            ),
        )
    end
    if 𝑁₁ > 0
        𝑁₁𝑁 = 𝑁₁ / 𝑁
        G₋ = -(geomean(abs.(a[a[:] .< 0, :])))
        varlogG₋ = var(log.(abs.(a[a[:] .< 0, :])))
        𝐸G₋ = abs(G₋) + (𝑁₁𝑁)^2 * varlogG₋ / (2 * 𝑁₁) * abs(G₋)
        t1 = (𝑁₁𝑁)^2 * (varlogG₋ / 𝑁₁)
        varG₋ = (𝑁₁𝑁)^2 * (varlogG₋ / 𝑁₁) * G₋
    else
        𝑁₁𝑁 = 0
        G₋ = 0
        varlogG₋ = 0
        𝐸G₋ = 0
        t1 = 0
        varG₋ = 0
    end
    if 𝑁₂ > 0
        𝑁₂𝑁 = 𝑁₂ / 𝑁
        G₊ = geomean(a[a[:] .> 0, :])
        varlogG₊ = var(log.(abs.(a[a[:] .> 0, :])))
        𝐸G₊ = G₊ + (𝑁₂𝑁)^2 * varlogG₊ / (2 * 𝑁₂) * (G₊)^2
        t2 = (𝑁₂𝑁)^2 * (varlogG₊ / 𝑁₂)
        varG₊ = (𝑁₂𝑁)^2 * (varlogG₊ / 𝑁₂) * (G₊)^2
    else
        𝑁₂𝑁 = 0
        G₊ = 0
        varlogG₊ = 0
        𝐸G₊ = 0
        t2 = 0
        varG₊ = 0
    end
    G = (𝑁₁ * G₋ + 𝑁₂ * G₊) / 𝑁
    𝐸G = G + (G / 2) * (t1 + t2)
    if 𝑁₁ == 0 || 𝑁₂ == 0
        covG = 0
    else
        covG = 𝐸G - 𝐸G₋ * 𝐸G₊
    end
    return (𝑁₁𝑁)^2 * varG₋ + (𝑁₂𝑁)^2 * varG₊ - ((2 * 𝑁₁ * 𝑁₂) / 𝑁^2) * covG
end
