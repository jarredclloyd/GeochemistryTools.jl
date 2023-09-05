#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-08-22
Edited: 2023-08-22

This source file contains functions to compute a line of best fit and its standard errors using the 'Mahon'
errors-in-variables regression algorithm of Mahon (1996) corrected by Stephan and Trappitsch (2023).

=#
# function exports
export eivlr_mahon

# caller functions
"""
    _eivlr_mahon(df::AbstractDataFrame; [se_level_in::Int=2, se_level_out::Int=2, se_type::AbstractString="abs",
        initial::Any=nothing])

Compute line of best fit using the "Mahon" errors-in-variables linear regression algorithm.

Input df as a DataFrame of 4 of 5 columns wide with column order (X, sX, Y, sY, [ρXY]).

# Keywords

  - `se_level_in::Int`: Standard error level of input data. Provide as a positive integer.
  - `se_level_out::Int`: Standard error level of output data. Provide as a positive integer.
  - `se_type::AbstractString`: Standard error type as a string of value `"abs"` OR `"rel"`. Values equal to
    `'a'`, `"absolute"`, `'r'`, and `"relative"` will also work. Case insensitive.
  - `initial::Any`: A value for the y-intercept. Can be input as a string key from an appropriate dictionary, as a single
    numeric value, or as a vector of the initial and its standard error (same `se_level_in` as input data). E.g. initial =
    "MDCInv", initial = 0.72, OR initial = [0.72, 0.01].

      + Dictionaries available are `dict_sr87_sr86i`
      + For a full list of available keys ∈ any dictionary type `keys(<dict_name>)`

# References

Mahon, KI (1996) 'The New “York” Regression: Application of an Improved Statistical Method to Geochemistry',
*International Geology Review*, 38(4):293–303. https://doi.org/10.1080/00206819709465336

Stephan, T & Trappitsch, R (2023) 'Reliable uncertainties: Error correlation, rotated error bars, and linear
regressions ∈ three-isotope plots and beyond', *International Journal of Mass Spectrometry*, 491:117053.
https://doi.org/10.1016/j.ijms.2023]].117053
"""
function _eivlr_mahon(X::AbstractArray, sX::AbstractArray, Y::AbstractArray, sY::AbstractArray, ρXY = nothing)
    nX::Int = length(X)
    if ρXY === nothing
        ρXY::AbstractArray = zeros(nX)
    elseif length(ρXY) !== nX
        ρXY = push!(zeros(nX - length(ρXY)))
    end
    β₀::AbstractFloat, β₁::AbstractFloat = coeffs(Polynomials.fit(X, Y, 1))
    βₑ::AbstractFloat = β₁
    σxy::AbstractArray = ρXY .* sX .* sY
    Ω::AbstractArray = 1 ./ (sY .^ 2 .+ β₁ .^ 2 .* sX .^ 2 .- 2 .* β₁ .* σxy)
    X̄::AbstractFloat = sum(Ω .* X) / sum(Ω)
    Ȳ::AbstractFloat = sum(Ω .* Y) / sum(Ω)
    U::AbstractArray = X .- X̄
    V::AbstractArray = Y .- Ȳ
    β₁ =
        sum(Ω .^ 2 .* V .* (U .* sY .^ 2 .+ β₁ .* V .* sX .^ 2 .- V .* σxy)) /
        sum(Ω .^ 2 .* U .* (U .* sY .^ 2 .+ β₁ .* V .* sX .^ 2 .- β₁ .* U .* σxy))
    n_iterations::Int = 1
    while abs(βₑ - β₁) > eps() && n_iterations < 1e6
        βₑ = β₁
        β₁ =
            sum(Ω .^ 2 .* V .* (U .* sY .^ 2 .+ β₁ .* V .* sX .^ 2 .- V .* σxy)) /
            sum(Ω .^ 2 .* U .* (U .* sY .^ 2 .+ β₁ .* V .* sX .^ 2 .- β₁ .* U .* σxy))
        n_iterations += 1
    end
    β₀ = Ȳ - β₁ * X̄
    # derivative calculations
    δθδβ₁ = _δθδβ₁(β₁, Ω, U, V, sX, sY, σxy)
    δθδX = zeros(AbstractFloat, length(X))
    Threads.@threads for i ∈ eachindex(X)
        δθδX[i] = _δθδXᵢ(i, β₁, Ω, U, V, sX, sY, σxy)
    end
    δθδY = zeros(AbstractFloat, length(X))
    Threads.@threads for i ∈ eachindex(X)
        δθδY[i] = _δθδYᵢ(i, β₁, Ω, U, V, sX, sY, σxy)
    end
    δβ₀δX = zeros(AbstractFloat, length(X))
    Threads.@threads for i ∈ eachindex(X)
        δβ₀δX[i] = _δβ₀δXᵢ(i, X̄, β₁, Ω, δθδX, δθδβ₁)
    end
    δβ₀δY = zeros(AbstractFloat, length(X))
    Threads.@threads for i ∈ eachindex(X)
        δβ₀δY[i] = _δβ₀δYᵢ(i, X̄, Ω, δθδY, δθδβ₁)
    end
    # variance calculations
    σβ₁² = sum(δθδX .^ 2 .* sX .^ 2 .+ δθδY .^ 2 .* sY .^ 2 .+ 2 .* σxy .* δθδX .* δθδY) / δθδβ₁^2
    σβ₀² = sum(δβ₀δX .^ 2 .* sX .^ 2 .+ δβ₀δY .^ 2 .* sY .^ 2 .+ 2 .* σxy .* δβ₀δX .* δβ₀δY)
    # X_intercept and variance {β₁ → 1/β₁; Xᵢ ⇄ Yᵢ; sX ⇄ sY; X̄ ⇄ Ȳ}
    # if calc_X_intercept != false
    X_intercept = -β₀ / β₁
    Ω = 1 ./ (sX .^ 2 .+ (1 / β₁) .^ 2 .* sY .^ 2 .- 2 .* (1 / β₁) .* σxy)
    δθδβ₁ = _δθδβ₁((1 / β₁), Ω, V, U, sY, sX, σxy)
    δθδX = zeros(AbstractFloat, length(X))
    Threads.@threads for i ∈ eachindex(X)
        δθδX[i] = _δθδXᵢ(i, (1 / β₁), Ω, V, U, sY, sX, σxy)
    end
    δθδY = zeros(AbstractFloat, length(X))
    Threads.@threads for i ∈ eachindex(X)
        δθδY[i] = _δθδYᵢ(i, (1 / β₁), Ω, V, U, sY, sX, σxy)
    end
    δβ₀δX = zeros(AbstractFloat, length(X))
    Threads.@threads for i ∈ eachindex(X)
        δβ₀δX[i] = _δβ₀δXᵢ(i, Ȳ, (1 / β₁), Ω, δθδX, δθδβ₁)
    end
    δβ₀δY = zeros(AbstractFloat, length(X))
    Threads.@threads for i ∈ eachindex(X)
        δβ₀δY[i] = _δβ₀δYᵢ(i, Ȳ, Ω, δθδY, δθδβ₁)
    end
    σX_intercept² = sum(δβ₀δX .^ 2 .* sX .^ 2 .+ δβ₀δY .^ 2 .* sY .^ 2 .+ 2 .* σxy .* δβ₀δX .* δβ₀δY)
    # end
    β₀SE = √(σβ₀²)
    β₁SE = √(σβ₁²)
    X_interceptSE = √(σX_intercept²)
    χ²::AbstractFloat = sum(Ω .* (Y .- β₁ .* X .- β₀) .^ 2)
    ν::Int = nX - 2
    χ²ᵣ::AbstractFloat = χ² / ν
    pval::AbstractFloat = ccdf(Chisq(ν), χ²)
    return Mahon(β₀, β₀SE, β₁, β₁SE, X_intercept, X_interceptSE, χ²ᵣ, pval, nX)
end

function _eivlr_mahon_fixedpoint(
    X::AbstractArray,
    sX::AbstractArray,
    Y::AbstractArray,
    sY::AbstractArray,
    ρXY = nothing;
    X₀Y₀::Tuple{Real,Real,Real,Real} = (0, 0, 0, 0),
)
    nX::Int = length(X)
    if ρXY === nothing
        ρXY::AbstractArray = zeros(nX)
    elseif length(ρXY) !== nX
        ρXY = push!(zeros(nX - length(ρXY)))
    end
    β₀::AbstractFloat, β₁::AbstractFloat = coeffs(Polynomials.fit(X, Y, 1))
    βₑ::AbstractFloat = β₁
    σxy::AbstractArray = ρXY .* sX .* sY
    Ω::AbstractArray = 1 ./ (sY .^ 2 .+ β₁ .^ 2 .* sX .^ 2 .- 2 .* β₁ .* σxy)
    X̄::AbstractFloat = sum(Ω .* X) / sum(Ω)
    Ȳ::AbstractFloat = sum(Ω .* Y) / sum(Ω)
    U::AbstractArray = X .- X̄
    V::AbstractArray = Y .- Ȳ
    β₁ =
        sum(Ω .^ 2 .* V .* (U .* sY .^ 2 .+ β₁ .* V .* sX .^ 2 .- V .* σxy)) /
        sum(Ω .^ 2 .* U .* (U .* sY .^ 2 .+ β₁ .* V .* sX .^ 2 .- β₁ .* U .* σxy))
    n_iterations::Int = 1
    while abs(βₑ - β₁) > eps() && n_iterations < 1e6
        βₑ = β₁
        β₁ =
            sum(Ω .^ 2 .* V .* (U .* sY .^ 2 .+ β₁ .* V .* sX .^ 2 .- V .* σxy)) /
            sum(Ω .^ 2 .* U .* (U .* sY .^ 2 .+ β₁ .* V .* sX .^ 2 .- β₁ .* U .* σxy))
        n_iterations += 1
    end
    β₀ = Ȳ - β₁ * X̄
    # derivative calculations
    δθδβ₁ = _δθδβ₁_fp(β₁, Ω, U, V, sX, sY, σxy)
    δθδX = _δθδXᵢ_fp(β₁, Ω, U, V, sX, sY, σxy)
    δθδY = _δθδYᵢ_fp(β₁, Ω, U, V, sX, sY, σxy)
    δβ₀δX = _δβ₀δXᵢ_fp(X̄, δθδX, δθδβ₁)
    δβ₀δY = _δβ₀δYᵢ_fp(X̄, δθδY, δθδβ₁)
    # variance calculations
    σβ₁² = sum(δθδX .^ 2 .* sX .^ 2 .+ δθδY .^ 2 .* sY .^ 2 .+ 2 .* σxy .* δθδX .* δθδY) / δθδβ₁^2
    σβ₀² = σβ₁² * X₀Y₀[1]
    # X_intercept and variance {β₁ → 1/β₁; Xᵢ ⇄ Yᵢ; sX ⇄ sY; X̄ ⇄ Ȳ}
    # if calc_X_intercept != false
    X_intercept = -β₀ / β₁
    Ω = 1 ./ (sX .^ 2 .+ (1 / β₁) .^ 2 .* sY .^ 2 .- 2 .* (1 / β₁) .* σxy)
    δθδβ₁ = _δθδβ₁_fp((1 / β₁), Ω, V, U, sY, sX, σxy)
    δθδX = _δθδXᵢ_fp((1 / β₁), Ω, V, U, sY, sX, σxy)
    δθδY = _δθδYᵢ_fp((1 / β₁), Ω, V, U, sY, sX, σxy)
    δβ₀δX = _δβ₀δXᵢ_fp(Ȳ, δθδX, δθδβ₁)
    δβ₀δY = _δβ₀δYᵢ_fp(Ȳ, δθδY, δθδβ₁)
    σX_intercept² = sum(δβ₀δX .^ 2 .* sX .^ 2 .+ δβ₀δY .^ 2 .* sY .^ 2 .+ 2 .* σxy .* δβ₀δX .* δβ₀δY)
    # end
    β₀SE = √(σβ₀²)
    β₁SE = √(σβ₁²)
    X_interceptSE = √(σX_intercept²)
    χ²::AbstractFloat = sum(Ω .* (Y .- β₁ .* X .- β₀) .^ 2)
    ν::Int = nX - 1
    χ²ᵣ::AbstractFloat = χ² / ν
    pval::AbstractFloat = ccdf(Chisq(ν), χ²)
    return Mahon(β₀, β₀SE, β₁, β₁SE, X_intercept, X_interceptSE, χ²ᵣ, pval, nX)
end

function _kroneckerδ(i::Integer, j::Integer)
    return i == j ? 1 : 0
end

function _δθδβ₁(
    β₁::AbstractFloat,
    Ω::AbstractArray,
    U::AbstractArray,
    V::AbstractArray,
    sX::AbstractArray,
    sY::AbstractArray,
    σxy::AbstractArray,
)
    return sum(Ω .^ 2 .* (2β₁ .* (U .* V .* sX .^ 2 .- U .^ 2 .* σxy) .+ (U .^ 2 .* sY .^ 2 .- V .^ 2 .* sX .^ 2))) +
           4 * sum(
               Ω .^ 3 .* (σxy .- β₁ .* sX .^ 2) .* (
                   β₁^2 .* (U .* V .* sX .^ 2 .- U .^ 2 .* σxy) .+ β₁ .* (U .^ 2 .* sY .^ 2 .- V .^ 2 .* sX .^ 2) .-
                   (U .* V .* sY .^ 2 .- V .^ 2 .* σxy)
               ),
           ) +
           2 *
           sum(Ω .^ 2 .* (-(β₁^2) .* U .* sX .^ 2 .+ 2β₁ .* V .* sX .^ 2 .+ U .* sY .^ 2 .- 2 .* V .* σxy)) *
           sum(Ω .^ 2 .* V .* (σxy .- β₁ .* sX .^ 2)) / sum(Ω) +
           2 *
           sum(Ω .^ 2 .* (-(β₁^2) .* V .* sX .^ 2 .+ 2(β₁^2) .* U .* σxy .- 2β₁ .* U .* sY .^ 2 .+ V .* sY .^ 2)) *
           sum(Ω .^ 2 .* U .* (σxy .- β₁ .* sX .^ 2)) / sum(Ω)
end

function _δθδXᵢ(
    i::Integer,
    β₁::AbstractFloat,
    Ω::AbstractArray,
    U::AbstractArray,
    V::AbstractArray,
    sX::AbstractArray,
    sY::AbstractArray,
    σxy::AbstractArray,
)
    tempⱼ = zeros(AbstractFloat, length(Ω))
    @simd for j ∈ eachindex(Ω)
        tempⱼ[j] =
            Ω[j]^2 *
            (_kroneckerδ(i, j) - Ω[i] / sum(Ω)) *
            (β₁^2 * V[j] * sX[j]^2 - 2 * β₁^2 * U[j] * σxy[j] + 2 * β₁ * U[j] * sY[j]^2 - V[j] * sY[j]^2)
    end
    δθδXᵢ = sum(tempⱼ)
    return δθδXᵢ
end

function _δθδYᵢ(
    i::Integer,
    β₁::AbstractFloat,
    Ω::AbstractArray,
    U::AbstractArray,
    V::AbstractArray,
    sX::AbstractArray,
    sY::AbstractArray,
    σxy::AbstractArray,
)
    tempⱼ = zeros(AbstractFloat, length(Ω))
    @simd for j ∈ eachindex(Ω)
        tempⱼ[j] =
            Ω[j]^2 *
            (_kroneckerδ(i, j) - Ω[j] / sum(Ω)) *
            (β₁^2 * U[j] * sX[j]^2 - 2 * β₁ * V[j] * sX[j]^2 - U[j] * sY[j]^2 + 2 * V[j] * σxy[j])
    end
    δθδYᵢ = sum(tempⱼ)
    return δθδYᵢ
end

function _δβ₀δXᵢ(
    i::Integer,
    X̄::AbstractFloat,
    β₁::AbstractFloat,
    Ω::AbstractArray,
    δθδX::AbstractArray,
    δθδβ₁::AbstractFloat,
)
    return -β₁ * Ω[i] / sum(Ω) - X̄ * δθδX[i] / δθδβ₁
end

function _δβ₀δYᵢ(i::Integer, X̄::AbstractFloat, Ω::AbstractArray, δθδY::AbstractArray, δθδβ₁::AbstractFloat)
    return Ω[i] / sum(Ω) - X̄ * δθδY[i] / δθδβ₁
end

# fixed point Mahon
function _δθδβ₁_fp(
    β₁::AbstractFloat,
    Ω::AbstractArray,
    U::AbstractArray,
    V::AbstractArray,
    sX::AbstractArray,
    sY::AbstractArray,
    σxy::AbstractArray,
)
    return sum(Ω .^ 2 .* (2β₁ .* (U .* V .* sX .^ 2 .- U .^ 2 .* σxy) .+ (U .^ 2 .* sY .^ 2 .- V .^ 2 .* sX .^ 2))) +
           4 * sum(
        Ω .^ 3 .* (σxy .- β₁ .* sX .^ 2) .* (
            β₁^2 .* (U .* V .* sX .^ 2 .- U .^ 2 .* σxy) .+ β₁ .* (U .^ 2 .* sY .^ 2 .- V .^ 2 .* sX .^ 2) .-
            (U .* V .* sY .^ 2 .- V .^ 2 .* σxy)
        ),
    )
end

function _δθδXᵢ_fp(
    β₁::AbstractFloat,
    Ω::AbstractArray,
    U::AbstractArray,
    V::AbstractArray,
    sX::AbstractArray,
    sY::AbstractArray,
    σxy::AbstractArray,
)
    return Ω .^ 2 .* (β₁^2 .* V .* sX .^ 2 .- 2 .* β₁^2 .* U .* σxy .+ 2 .* β₁ .* U .* sY .^ 2 .- V .* sY .^ 2)
end

function _δθδYᵢ_fp(
    β₁::AbstractFloat,
    Ω::AbstractArray,
    U::AbstractArray,
    V::AbstractArray,
    sX::AbstractArray,
    sY::AbstractArray,
    σxy::AbstractArray,
)
    return Ω .^ 2 .* (β₁^2 .* U .* sX^0.2 .- 2 .* β₁ .* V .* sX .^ 2 .- U .* sY .^ 2 .+ 2 .* V .* σxy)
end

function _δβ₀δXᵢ_fp(X̄::AbstractFloat, δθδX::AbstractArray, δθδβ₁::AbstractFloat)
    return -X̄ * δθδX / δθδβ₁
end

function _δβ₀δYᵢ_fp(X̄::AbstractFloat, δθδY::AbstractArray, δθδβ₁::AbstractFloat)
    return -X̄ * δθδY / δθδβ₁
end
