#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-08-22
Edited: 2023-08-22

This source file contains functions to compute a line of best fit and its standard errors using the 'Mahon'
errors-in-variables regression algorithm of Mahon (1996) corrected by Stephan and Trappitsch (2023).

=#
# function exports
export fit_mahon

# caller functions
"""
    fit_mahon(df::AbstractDataFrame; [se_level_in::Int=2, se_level_out::Int=2, se_type::AbstractString="abs",
        initial::Any=nothing])

Compute line of best fit using the "Mahon" errors-in-variables linear regression algorithm.

Input df as a DataFrame of 4 of 5 columns wide with column order (x, σx, y, σy, [ρxy]).

# Keywords

  - `se_level_in::Int`: Standard error level of input data. Provide as a positive integer.

  - `se_level_out::Int`: Standard error level of output data. Provide as a positive integer.
  - `se_type::AbstractString`: Standard error type as a string of value `"abs"` OR `"rel"`. Values equal to
    `'a'`, `"absolute"`, `'r'`, and `"relative"` will also work. Case insensitive.
  - `initial::Any`: A value for the y-intercept. Can be input as a string key from an appropriate dictionary, as a single
    numeric value, or as a vector of the initial and its standard error (same `se_level_in` as input data). E.g. initial =
    "MDCInv", initial = 0.72, OR initial = [0.72, 0.01].

      + Dictionaries available are `INITIAL_Sr`
      + For a full list of available keys ∈ any dictionary type `keys(<dict_name>)`

# References

Mahon, KI (1996) 'The New “York” Regression: Application of an Improved Statistical Method to Geochemistry',
*International Geology Review*, 38(4):293–303. https://doi.org/10.1080/00206819709465336

Stephan, T & Trappitsch, R (2023) 'Reliable uncertainties: Error correlation, rotated error bars, and linear
regressions ∈ three-isotope plots and beyond', *International Journal of Mass Spectrometry*, 491:117053.
https://doi.org/10.1016/j.ijms.2023]].117053
"""
function fit_mahon(
    df::AbstractDataFrame;
    se_level_in::Int = 2,
    se_level_out::Int = 2,
    se_type::AbstractString = "abs",
    initial::Any = nothing,
) end

"""
    _eivlr_mahon(x::Vector{<:Real}, σx::Vector{<:Real}, y::Vector{<:Real}, σy::Vector{<:Real}, ρxy::Union{Nothing, Vector{<:Real}} = nothing)

Compute line of best fit using the "Mahon" errors-in-variables regression algorithm.
"""
function _eivlr_mahon(
    x::Vector{<:Real},
    σx::Vector{<:Real},
    y::Vector{<:Real},
    σy::Vector{<:Real},
    ρxy::Union{Nothing, Vector{<:Real}} = nothing,
)
    if _check_equal_length(x, σx, y, σy) != true
        throw(ArgumentError("The length of x, σx, y and σy must be the same."))
    end
    𝑁::Int = length(x)
    if ρxy === nothing
        ρxy::Vector{Float64} = zeros(𝑁)
    elseif length(ρxy) !== 𝑁
        throw(ArgumentError("The length of ρxy must be the same as x, σx, y and σy."))
    end

    # initial slope and intercept from (OLS)
    β₀::Float64, β₁::Float64 = hcat(ones(𝑁), x) \ y
    βₑ::Float64 = β₁

    # initial fit via Mahon method
    σxy::Vector{Float64} = ρxy .* σx .* σy
    Ω::Vector{Float64} = @. (1 / (σy^2 + β₁^2 * σx^2 - 2 * β₁ * σxy))
    x̄::Float64 = sum(Ω .* x) / sum(Ω)
    ȳ::Float64 = sum(Ω .* y) / sum(Ω)
    u::Vector{Float64} = x .- x̄
    v::Vector{Float64} = y .- ȳ
    β₁ =
        sum(@.(Ω^2 * v * (u * σy^2 + β₁ * v * σx^2 - v * σxy))) /
        sum(@.(Ω^2 * u * (u * σy^2 + β₁ * v * σx^2 - β₁ * u * σxy)))

    n_iterations::Integer = 1
    # iterative solve via Mahon method
    while abs(βₑ - β₁) > 1e-15 && n_iterations < 1e6
        βₑ = β₁
        Ω = @.(1 / (σy^2 + βₑ^2 * σx^2 - 2 * βₑ * σxy))
        x̄ = sum(Ω .* x) / sum(Ω)
        ȳ = sum(Ω .* y) / sum(Ω)
        u = x .- x̄
        v = y .- ȳ
        β₁ =
            sum(@. (Ω^2 * v * (u * σy^2 + β₁ * v * σx^2 - v * σxy))) /
            sum(@. (Ω^2 * u * (u * σy^2 + β₁ * v * σx^2 - β₁ * u * σxy)))
        n_iterations += 1
    end
    β₀ = ȳ - β₁ * x̄
    χ²::Float64 = sum(@. (Ω * (y - β₁ * x - β₀)^2))
    ν::Int = 𝑁 > 2 ? 𝑁 - 2 : 1
    χ²ᵣ::Float64 = χ² / ν
    pval::Float64 = ccdf(Chisq(ν), χ²)
    # derivative calculations
    δθδβ₁ = _δθδβ₁(β₁, Ω, u, v, σx, σy, σxy)
    δθδX = zeros(Float64, 𝑁)
    Threads.@threads for i ∈ eachindex(x)
        δθδX[i] = _δθδXᵢ(i, β₁, Ω, u, v, σx, σy, σxy)
    end
    δθδY = zeros(Float64, 𝑁)
    Threads.@threads for i ∈ eachindex(x)
        δθδY[i] = _δθδYᵢ(i, β₁, Ω, u, v, σx, σy, σxy)
    end
    δβ₀δX = zeros(Float64, 𝑁)
    Threads.@threads for i ∈ eachindex(x)
        δβ₀δX[i] = _δβ₀δXᵢ(i, x̄, β₁, Ω, δθδX, δθδβ₁)
    end
    δβ₀δY = zeros(Float64, 𝑁)
    Threads.@threads for i ∈ eachindex(x)
        δβ₀δY[i] = _δβ₀δYᵢ(i, x̄, Ω, δθδY, δθδβ₁)
    end
    # variance calculations
    σβ₁² = sum(@.(δθδX^2 * σx^2 + δθδY^2 * σy^2 + 2 * σxy * δθδX * δθδY)) / δθδβ₁^2
    σβ₀² = sum(@.(δβ₀δX^2 * σx^2 + δβ₀δY^2 * σy^2 + 2 * σxy * δβ₀δX * δβ₀δY))
    X_intercept = -β₀ / β₁
    #= X_intercept and variance {β₁ → 1/β₁; Xᵢ ⇄ Yᵢ; σx ⇄ σy; x̄ ⇄ ȳ}
    Ω = @.(1 / (σx^2 + (1 / β₁)^2 * σy^2 - 2 * (1 / β₁) * σxy))
    δθδβ₁ = _δθδβ₁((1 / β₁), Ω, u, v, σy, σx, σxy)
    δθδX = zeros(Float64, 𝑁)
    Threads.@threads for i ∈ eachindex(x)
        δθδX[i] = _δθδXᵢ(i, (1 / β₁), Ω, u, v, σy, σx, σxy)
    end
    δθδY = zeros(Float64, 𝑁)
    Threads.@threads for i ∈ eachindex(x)
        δθδY[i] = _δθδYᵢ(i, (1 / β₁), Ω, u, v, σy, σx, σxy)
    end
    δβ₀δX = zeros(Float64, 𝑁)
    Threads.@threads for i ∈ eachindex(x)
        δβ₀δX[i] = _δβ₀δXᵢ(i, ȳ, (1 / β₁), Ω, δθδX, δθδβ₁)
    end
    δβ₀δY = zeros(Float64, 𝑁)
    Threads.@threads for i ∈ eachindex(x)
        δβ₀δY[i] = _δβ₀δYᵢ(i, ȳ, Ω, δθδY, δθδβ₁)
    end
    σX_intercept² = sum(@.(δβ₀δX^2 * σy^2 + δβ₀δY^2 * σx^2 + 2 * σxy * δβ₀δX * δβ₀δY))
    =#
    β₀SE = √(σβ₀²)
    β₁SE = √(σβ₁²)
    σᵦ₁ᵦ₀::Float64 = -x̄ * σβ₁²
    X_interceptSE = sqrt((β₀SE / β₀)^2 + (β₁SE / β₁)^2 - 2 * σᵦ₁ᵦ₀ / (β₀ * β₁)) #=√(σX_intercept²)=#
    return MahonNonFixed(
        β₀,
        β₀SE,
        β₁,
        β₁SE,
        X_intercept,
        X_interceptSE,
        χ²ᵣ,
        pval,
        σᵦ₁ᵦ₀,
        𝑁
    )
end

function _eivlr_mahon_fixedpoint(
    x::AbstractArray,
    σx::AbstractArray,
    y::AbstractArray,
    σy::AbstractArray,
    ρxy = nothing;
    x₀y₀::Tuple{Real,Real,Real,Real} = (0, 0, 0, 0),
)
    if _check_equal_length(x, σx, y, σy) != true
        throw(ArgumentError("The length of x, σx, y and σy must be the same."))
    end
    𝑁::Int = length(x)
    if ρxy === nothing
        ρxy::Vector{Float64} = zeros(𝑁)
    elseif length(ρxy) !== 𝑁
        throw(ArgumentError("The length of ρxy must be the same as x, σx, y and σy."))
    end
    β₀::Float64, β₁::Float64 = hcat(ones(𝑁), x) \ y
    βₑ::Float64 = β₁
    σxy::Vector{Float64} = ρxy .* σx .* σy
    Ω::Vector{Float64} = @. (1 / (σy^2 + β₁^2 * σx^2 - 2 * β₁ * σxy))
    x̄::Float64 = sum(Ω .* x) / sum(Ω)
    ȳ::Float64 = sum(Ω .* y) / sum(Ω)
    u::Vector{Float64} = x .- x̄
    v::Vector{Float64} = y .- ȳ
    β₁ =
        sum(@.(Ω^2 * v * (u * σy^2 + β₁ * v * σx^2 - v * σxy))) /
        sum(@.(Ω^2 * u * (u * σy^2 + β₁ * v * σx^2 - β₁ * u * σxy)))
    n_iterations::Integer = 1
    while abs(βₑ - β₁) > 1e-15 && n_iterations < 1e6
        βₑ = β₁
        Ω = @.(1 / (σy^2 + βₑ^2 * σx^2 - 2 * βₑ * σxy))
        x̄ = sum(Ω .* x) / sum(Ω)
        ȳ = sum(Ω .* y) / sum(Ω)
        u = x .- x̄
        v = y .- ȳ
        β₁ =
            sum(@. (Ω^2 * v * (u * σy^2 + β₁ * v * σx^2 - v * σxy))) /
            sum(@. (Ω^2 * u * (u * σy^2 + β₁ * v * σx^2 - β₁ * u * σxy)))
        n_iterations += 1
    end
    β₀ = ȳ - β₁ * x̄
    χ²::Float64 = sum(@.(Ω * (y - β₁ * x - β₀)^2))
    ν::Int = 𝑁 > 2 ? 𝑁 - 2 : 1
    χ²ᵣ::Float64 = χ² / ν
    pval::Float64 = ccdf(Chisq(ν), χ²)
    # derivative calculations
    δθδβ₁ = _δθδβ₁_fp(β₁, Ω, u, v, σx, σy, σxy)
    δθδX = _δθδXᵢ_fp(β₁, Ω, u, v, σx, σy, σxy)
    δθδY = _δθδYᵢ_fp(β₁, Ω, u, v, σx, σy, σxy)
    δβ₀δX = _δβ₀δXᵢ_fp(x̄, δθδX, δθδβ₁)
    δβ₀δY = _δβ₀δYᵢ_fp(x̄, δθδY, δθδβ₁)
    # variance calculations
    σβ₁² = sum(@.(δθδX^2 * σx^2 + δθδY^2 * σy^2 + 2 * σxy * δθδX * δθδY)) / δθδβ₁^2
    σβ₀² = σβ₁² * x₀y₀[1]
    X_intercept = -β₀ / β₁
    #= X_intercept and variance {β₁ → 1/β₁; Xᵢ ⇄ Yᵢ; σx ⇄ σy; x̄ ⇄ ȳ}
    # if calc_X_intercept != false
    Ω = @.(1 / (σx^2 + (1 / β₁)^2 * σy^2 - 2 * (1 / β₁) * σxy))
    δθδβ₁ = _δθδβ₁_fp((1 / β₁), Ω, v, u, σy, σx, σxy)
    δθδX = _δθδXᵢ_fp((1 / β₁), Ω, v, u, σy, σx, σxy)
    δθδY = _δθδYᵢ_fp((1 / β₁), Ω, v, u, σy, σx, σxy)
    δβ₀δX = _δβ₀δXᵢ_fp(ȳ, δθδX, δθδβ₁)
    δβ₀δY = _δβ₀δYᵢ_fp(ȳ, δθδY, δθδβ₁)
    σX_intercept² = sum(@.(δβ₀δX^2 * σx^2 + δβ₀δY^2 * σy^2 + 2 * σxy * δβ₀δX * δβ₀δY))
    =#
    β₀SE = √(σβ₀²)
    β₁SE = √(σβ₁²)
    σᵦ₁ᵦ₀::Float64 = -x̄ * σβ₁²
    X_interceptSE = sqrt((β₀SE / β₀)^2 + (β₁SE / β₁)^2 - 2 * σᵦ₁ᵦ₀ / (β₀ * β₁)) #=√(σX_intercept²)=#
    return MahonFixed(β₀, β₀SE, β₁, β₁SE, X_intercept, X_interceptSE, χ²ᵣ, pval, σᵦ₁ᵦ₀, 𝑁, x₀y₀)
end

function _kroneckerδ(i::Integer, j::Integer)
    return i == j ? 1 : 0
end

function _δθδβ₁(
    β₁::AbstractFloat,
    Ω::AbstractArray,
    u::AbstractArray,
    v::AbstractArray,
    σx::AbstractArray,
    σy::AbstractArray,
    σxy::AbstractArray,
)
    return sum(@.(Ω^2 * (2β₁ * (u * v * σx^2 - u^2 * σxy) + (u^2 * σy^2 - v^2 * σx^2)))) +
           4 * sum(
               @.(
                   Ω^3 *
                   (σxy - β₁ * σx^2) *
                   (
                       β₁^2 * (u * v * σx^2 - u^2 * σxy) + β₁ * (u^2 * σy^2 - v^2 * σx^2) -
                       (u * v * σy^2 - v^2 * σxy)
                   )
               ),
           ) +
           2 *
           sum(@.(Ω^2 * (-(β₁^2) * u * σx^2 + 2β₁ * v * σx^2 + u * σy^2 - 2 * v * σxy))) *
           sum(@.(Ω^2 * v * (σxy - β₁ * σx^2))) / sum(Ω) +
           2 *
           sum(
               @.(
                   Ω^2 *
                   (-(β₁^2) * v * σx^2 + 2(β₁^2) * u * σxy - 2β₁ * u * σy^2 + v * σy^2),
               )
           ) *
           sum(@.(Ω^2 * u * (σxy - β₁ * σx^2))) / sum(Ω)
end

function _δθδXᵢ(
    i::Integer,
    β₁::AbstractFloat,
    Ω::AbstractArray,
    u::AbstractArray,
    v::AbstractArray,
    σx::AbstractArray,
    σy::AbstractArray,
    σxy::AbstractArray,
)
    δθδXᵢ = zeros(Float64, length(Ω))
    @simd for j ∈ eachindex(Ω)
        δθδXᵢ[j] =
            Ω[j]^2 *
            (_kroneckerδ(i, j) - Ω[i] / sum(Ω)) *
            (
                β₁^2 * v[j] * σx[j]^2 - 2 * β₁^2 * u[j] * σxy[j] + 2 * β₁ * u[j] * σy[j]^2 -
                v[j] * σy[j]^2
            )
    end
    δθδXᵢ = sum(δθδXᵢ)
    return δθδXᵢ
end

function _δθδYᵢ(
    i::Integer,
    β₁::AbstractFloat,
    Ω::AbstractArray,
    u::AbstractArray,
    v::AbstractArray,
    σx::AbstractArray,
    σy::AbstractArray,
    σxy::AbstractArray,
)
    δθδYᵢ = zeros(Float64, length(Ω))
    @simd for j ∈ eachindex(Ω)
        δθδYᵢ[j] =
            Ω[j]^2 *
            (_kroneckerδ(i, j) - Ω[i] / sum(Ω)) *
            (
                β₁^2 * u[j] * σx[j]^2 - 2 * β₁ * v[j] * σx[j]^2 - u[j] * σy[j]^2 +
                2 * v[j] * σxy[j]
            )
    end
    δθδYᵢ = sum(δθδYᵢ)
    return δθδYᵢ
end

function _δβ₀δXᵢ(
    i::Integer,
    x̄::AbstractFloat,
    β₁::AbstractFloat,
    Ω::AbstractArray,
    δθδX::AbstractArray,
    δθδβ₁::AbstractFloat,
)
    return -β₁ * Ω[i] / sum(Ω) - x̄ * δθδX[i] / δθδβ₁
end

function _δβ₀δYᵢ(
    i::Integer,
    x̄::AbstractFloat,
    Ω::AbstractArray,
    δθδY::AbstractArray,
    δθδβ₁::AbstractFloat,
)
    return Ω[i] / sum(Ω) - x̄ * δθδY[i] / δθδβ₁
end

# fixed point Mahon
function _δθδβ₁_fp(
    β₁::AbstractFloat,
    Ω::AbstractArray,
    u::AbstractArray,
    v::AbstractArray,
    σx::AbstractArray,
    σy::AbstractArray,
    σxy::AbstractArray,
)
    return sum(@.(Ω^2 * (2β₁ * (u * v * σx^2 - u^2 * σxy) + (u^2 * σy^2 - v^2 * σx^2)))) +
           4 * sum(
        @.(
            Ω^3 *
            (σxy - β₁ * σx^2) *
            (
                β₁^2 * (u * v * σx^2 - u^2 * σxy) + β₁ * (u^2 * σy^2 - v^2 * σx^2) -
                (u * v * σy^2 - v^2 * σxy)
            ),
        )
    )
end

function _δθδXᵢ_fp(
    β₁::AbstractFloat,
    Ω::AbstractArray,
    u::AbstractArray,
    v::AbstractArray,
    σx::AbstractArray,
    σy::AbstractArray,
    σxy::AbstractArray,
)
    return @.(Ω^2 * (β₁^2 * v * σx^2 - 2 * β₁^2 * u * σxy + 2 * β₁ * u * σy^2 - v * σy^2))
end

function _δθδYᵢ_fp(
    β₁::AbstractFloat,
    Ω::AbstractArray,
    u::AbstractArray,
    v::AbstractArray,
    σx::AbstractArray,
    σy::AbstractArray,
    σxy::AbstractArray,
)
    return @.(Ω^2 * (β₁^2 * u * σx^2 - 2 * β₁ * v * σx^2 - u * σy^2 + 2 * v * σxy))
end

function _δβ₀δXᵢ_fp(x̄::AbstractFloat, δθδX::AbstractArray, δθδβ₁::AbstractFloat)
    return -x̄ * δθδX / δθδβ₁
end

function _δβ₀δYᵢ_fp(x̄::AbstractFloat, δθδY::AbstractArray, δθδβ₁::AbstractFloat)
    return -x̄ * δθδY / δθδβ₁
end
