#= Preamble
This file contains functions for errors-in-variables linear regression.
=#

export york, yorkfit, deming, demingfit, jackknife, mahon

#Caller functions
"""
    demingfit(df::AbstractDataFrame)

Compute line of best fit using Deming errors-in-variables regression algorithm.

Input df as a DataFrame of 2 columns wide with order (X, Y). Standard errors for β₀ and β₁.
"""
function demingfit(df::AbstractDataFrame)
    β₀::AbstractFloat, β₁::AbstractFloat = deming(df[!, 1]::AbstractArray, df[!, 2]::AbstractArray)
    if isnan(β₀) || isnan(β₁)
        println("Solution not computed!")
    else
        β₀SE::AbstractFloat, β₁SE::AbstractFloat, n = jackknife(df, "deming")
    end
    return β₀, β₀SE, β₁, β₁SE, n
end

"""
    yorkfit(df::AbstractDataFrame; [se_level_in::Int=2, se_level_out::Int=2, se_type::AbstractString="abs",
        initial::Any=nothing])

Compute line of best fit using York errors-in-variables linear regression algorithm.

Input df as a DataFrame of 4 of 5 columns wide with column order (X, sX, Y, sY, [ρXY]).

# Keywords
- `se_level_in::Int`: Standard error level of input data. Provide as a positive integer.
- `se_level_out::Int`: Standard error level of output data. Provide as a positive integer.
- `se_type::AbstractString`: Standard error type as a string of value `"abs"` OR `"rel"`. Values equal to
    `'a'`, `"absolute"`, `'r'`, and `"relative"` will also work. Case insensitive.
- `initial::Any`: A value for the y-intercept. Can be input as a string key from an appropriate dictionary, as a single
    numeric value, or as a vector of the initial and its standard error (same `se_level_in` as input data). E.g. initial =
        "MDCInv", initial = 0.72, OR initial = [0.72, 0.01].
    - Dictionaries available are `dict_sr87_sr86i`
    - For a full list of available keys in any dictionary type `keys(<dict_name>)`

# Example
```julia-repl
julia> york(df, 2, "abs"; initial = "MDCInv")

```

# Reference
York D. et al. 2004 "Unified equations for the slope, intercept, and standard errors of the best straight line", *American
Journal of Physics*, 72(3), doi:https://doi.org/10.1119/1.1632486.
"""
function yorkfit(
    df::AbstractDataFrame;
    se_level_in::Int = 2,
    se_level_out::Int = 2,
    se_type::AbstractString = "abs",
    initial::Any = nothing,
)
    dfCols::Int = ncol(df)
    dfRows::Int = nrow(df)
    if initial !== nothing
        if isa(initial, String) == true && haskey(dict_sr87_sr86i, initial) == true
            initial = deepcopy(dict_sr87_sr86i[initial])
        elseif length(initial) == 1 && isa(initial, AbstractFloat)
            initial = [1e-10, 1e-18, initial, 0.01]
        elseif length(initial) == 2 && isa(initial[1], AbstractFloat) && isa(initial[2], AbstractFloat)
            initial = [1e-10, 1e-18, initial[1], initial[2]]
        elseif length(initial) >= 3
            throw(ArgumentError("Only the initial ratio value and its uncertainty is required."))
        else
            throw(ArgumentError("""
            initial variable is malformed.
            Use either a string that matches a fixed value or reference material name:
                e.g. "MDC", "MDCInv", or "SolarSystem"
            Or input a single numeric value or two-length numeric vector:
                e.g. 0.72 or [0.7200, 0.0010]

            For a full list of available reference material keys in a relevant dictionary: keys(<dict_name>)

            Dictionaries currently available are `dict_sr87_sr86i`
            """))
        end
        if dfCols == 4
            push!(df, initial)
        elseif dfCols == 5
            push!(initial, 0)
            push!(df, initial)
        end
    end
    if occursin('a', lowercase.(se_type)) == true ||
        occursin("abs", lowercase.(se_type)) == true ||
        occursin("absolute", lowercase.(se_type)) == true
        if dfCols == 5
            β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX = york(
                df[!, 1], df[!, 2] ./ se_level_in, df[!, 3], df[!, 4] ./ se_level_in, df[!, 5]
            )
        elseif dfCols == 4
            β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX = york(
                df[!, 1], df[!, 2] ./ se_level_in, df[!, 3], df[!, 4] ./ se_level_in
            )
        else
            throw(ArgumentError("Column width is not equal to 4 or 5. Some data is missing."))
        end
    elseif occursin('r', lowercase.(se_type)) == true ||
        occursin("rel", lowercase.(se_type)) == true ||
        occursin("relative", lowercase.(se_type)) == true
        if dfCols == 5
            β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX = york(
                df[!, 1],
                (df[!, 2] .* df[!, 1]) ./ se_level_in,
                df[!, 3],
                (df[!, 4] .* df[!, 3]) ./ se_level_in,
                df[!, 5],
            )
        elseif dfCols == 4
            β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX = york(
                df[!, 1], (df[!, 2] .* df[!, 1]) ./ se_level_in, df[!, 3], (df[!, 4] .* df[!, 3]) ./ se_level_in
            )
        else
            throw(ArgumentError("Column width is not equal to 4 or 5. Some data is missing."))
        end
    end
    if initial !== nothing
        popat!(df, dfRows + 1)
    end
    β₀SE *= se_level_out
    β₁SE *= se_level_out
    return β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX
end

#Base functions
function deming(X::AbstractArray, Y::AbstractArray)
    nX::Int = length(X)
    if nX !== length(Y)
        throw(
            ArgumentError("Lengths of data vectors are not equal. Please ensure all data (X & Y) are of equal length.")
        )
    else
        λ::AbstractFloat = (var(Y)) / (var(X))
        X̄::AbstractFloat = mean(X)
        Ȳ::AbstractFloat = mean(Y)
        Sxx::AbstractFloat = sum((X .- X̄) .^ 2)
        Syy::AbstractFloat = sum((Y .- Ȳ) .^ 2)
        Sxy::AbstractFloat = sum((X .- X̄) .* (Y .- Ȳ))
        β₁::AbstractFloat = (Syy - λ * Sxx + √((Syy - λ * Sxx)^2 + 4 * λ * Sxy^2)) / (2 * Sxy)
        β₀::AbstractFloat = Ȳ - β₁ * X̄
    end
    return β₀, β₁
end

#jackknife deming
function jackknife(df::AbstractDataFrame, method::AbstractString = "deming")
    nX::Int = nrow(df)
    if method == "deming"
        jackknifed_estimates::AbstractArray = zeros(AbstractFloat, nX, 2)
        @threads for j_position in 1:nX
            β₀est::AbstractFloat, β₁est::AbstractFloat = deming(df[Not(j_position), 1], df[Not(j_position), 2])
            jackknifed_estimates[j_position, 1] = β₀est
            jackknifed_estimates[j_position, 2] = β₁est
        end
        β₀SE::AbstractFloat = sem(jackknifed_estimates[:, 1])
        β₁SE::AbstractFloat = sem(jackknifed_estimates[:, 2])
        return β₀SE, β₁SE, nX
    end
end

function york(X::AbstractArray, sX::AbstractArray, Y::AbstractArray, sY::AbstractArray, ρXY = nothing)
    nX::Int = length(X)
    if ρXY === nothing
        ρXY::AbstractArray{AbstractFloat} = zeros(nX)
    elseif length(ρXY) !== nX
        ρXY = push!(zeros(nX - length(ρXY)))
    end
    β₀::AbstractFloat, β₁::AbstractFloat = coeffs(Polynomials.fit(X, Y, 1))
    βₑ::AbstractFloat = β₁
    ωXᵢ::AbstractArray{AbstractFloat} = 1 ./ sX .^ 2
    ωYᵢ::AbstractArray{AbstractFloat} = 1 ./ sY .^ 2
    α::AbstractArray{AbstractFloat} = .√(ωXᵢ .* ωYᵢ)
    Ω::AbstractArray{AbstractFloat} = ωXᵢ .* ωYᵢ ./ (ωXᵢ .+ β₁^2 .* ωYᵢ .- 2 .* β₁ .* ρXY .* α)
    X̄::AbstractFloat = sum(Ω .* X) / sum(Ω)
    Ȳ::AbstractFloat = sum(Ω .* Y) / sum(Ω)
    U::AbstractArray{AbstractFloat} = X .- X̄
    V::AbstractArray{AbstractFloat} = Y .- Ȳ
    βᵢ::AbstractArray{AbstractFloat} = Ω .* (U ./ ωYᵢ + β₁ * V ./ ωXᵢ - (β₁ * U + V) .* ρXY ./ α)
    β₁ = sum(Ω .* βᵢ .* V) / sum(Ω .* βᵢ .* U)
    n_iterations::Int = 1
    while (βₑ / β₁ - 1)^2 > eps() && n_iterations < 1e6
        Ω = ωXᵢ .* ωYᵢ ./ (ωXᵢ .+ β₁^2 .* ωYᵢ .- 2 .* β₁ .* ρXY .* α)
        X̄ = sum(Ω .* X) / sum(Ω)
        Ȳ = sum(Ω .* Y) / sum(Ω)
        U = X .- X̄
        V = Y .- Ȳ
        βᵢ = Ω .* (U ./ ωYᵢ + β₁ * V ./ ωXᵢ - (β₁ * U + V) .* ρXY ./ α)
        βₑ = β₁
        β₁ = sum(Ω .* βᵢ .* V) / sum(Ω .* βᵢ .* U)
        n_iterations = n_iterations + 1
    end
    β₀ = Ȳ - β₁ * X̄
    xᵢ = X̄ .+ βᵢ
    x̄ = sum(Ω .* xᵢ) / sum(Ω)
    υ::AbstractArray{AbstractFloat} = xᵢ .- x̄
    β₁SE::AbstractFloat = √(1 / sum(Ω .* υ .^ 2))
    β₀SE::AbstractFloat = √(1 / sum(Ω) + (x̄ * β₁SE)^2)
    σᵦ₁ᵦ₀::AbstractFloat = - x̄ * β₁SE^2
    χ²::AbstractFloat = sum(Ω .* (Y .- β₁ .* X .- β₀) .^ 2)
    ν::Int = nX - 2
    χ²ᵣ::AbstractFloat = χ² / ν
    pval::AbstractFloat = ccdf(Chisq(ν), χ²)
    return β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX, X̄, Ȳ
end

#=
function MonteCarloEIV()

end
=#

#= Mahon regression and associated functions
Mahon (1996) "new York regression" corrected by Stephan & Trappitsch (2023)
=#

function mahon(X::AbstractArray, sX::AbstractArray, Y::AbstractArray, sY::AbstractArray, ρXY = nothing)
    nX::Int = length(X)
    if ρXY === nothing
        ρXY::AbstractArray = zeros(nX)
    elseif length(ρXY) !== nX
        ρXY = push!(zeros(nX - length(ρXY)))
    end
    β₀::AbstractFloat, β₁::AbstractFloat = coeffs(Polynomials.fit(X, Y, 1))
    βₑ::AbstractFloat = β₁
    σxy::AbstractArray = ρXY .* sX .* sY
    Ω::AbstractArray = 1 ./ (sY .^2 .+ β₁ .^2 .* sX .^2 .- 2 .* β₁ .* σxy)
    X̄::AbstractFloat = sum(Ω .* X) / sum(Ω)
    Ȳ::AbstractFloat = sum(Ω .* Y) / sum(Ω)
    U::AbstractArray = X .- X̄
    V::AbstractArray = Y .- Ȳ
    β₁ =
        sum(Ω .^2 .* V .* (U .* sY .^2 .+ β₁ .* V .* sX .^2 .- V .* σxy)) /
        sum(Ω .^2 .* U .* (U .* sY .^2 .+ β₁ .* V .* sX .^2 .- β₁ .* U .* σxy))
    n_iterations::Int = 1
    while abs(βₑ - β₁) > eps() && n_iterations < 1e6
        βₑ = β₁
        β₁ =
            sum(Ω .^2 .* V .* (U .* sY .^2 .+ β₁ .* V .* sX .^2 .- V .* σxy)) /
            sum(Ω .^2 .* U .* (U .* sY .^2 .+ β₁ .* V .* sX .^2 .- β₁ .* U .* σxy))

        n_iterations += 1
    end
    β₀ = Ȳ - β₁ * X̄
    # derivative calculations
    δθδβ₁ = _δθδβ₁(β₁, Ω, U, V, sX, sY, σxy)
    δθδX = zeros(AbstractFloat, length(X))
    Threads.@threads for i in eachindex(X)
        δθδX[i] = _δθδXᵢ(i, β₁, Ω, U, V, sX, sY, σxy)
    end
    δθδY = zeros(AbstractFloat, length(X))
    Threads.@threads for i in eachindex(X)
        δθδY[i] = _δθδYᵢ(i, β₁, Ω, U, V, sX, sY, σxy)
    end
    δβ₀δX = zeros(AbstractFloat, length(X))
    Threads.@threads for i in eachindex(X)
        δβ₀δX[i] = _δβ₀δXᵢ(i, X̄, β₁, Ω, δθδX, δθδβ₁)
    end
    δβ₀δY = zeros(AbstractFloat, length(X))
    Threads.@threads for i in eachindex(X)
        δβ₀δY[i] = _δβ₀δYᵢ(i, X̄, Ω, δθδY, δθδβ₁)
    end
    # variance calculations
    σβ₁² = sum(δθδX .^2 .* sX .^2 .+ δθδY .^2 .* sY .^2 .+ 2 .* σxy .* δθδX .* δθδY) / δθδβ₁^2
    σβ₀² = sum(δβ₀δX .^ 2 .* sX .^ 2 .+ δβ₀δY .^ 2 .* sY .^ 2 .+ 2 .* σxy .* δβ₀δX .* δβ₀δY)
    # X_intercept and variance {β₁ → 1/β₁; Xᵢ ⇄ Yᵢ; sX ⇄ sY; X̄ ⇄ Ȳ}
    # if calc_X_intercept != false
        X_intercept = -β₀ / β₁
        Ω = 1 ./ (sX .^2 .+ (1/β₁) .^2 .* sY .^2 .- 2 .* (1/β₁) .* σxy)
        δθδβ₁ = _δθδβ₁((1/β₁), Ω, V, U, sY, sX, σxy)
        δθδX = zeros(AbstractFloat, length(X))
        Threads.@threads for i in eachindex(X)
            δθδX[i] = _δθδXᵢ(i, (1/β₁), Ω, V, U, sY, sX, σxy)
        end
        δθδY = zeros(AbstractFloat, length(X))
        Threads.@threads for i in eachindex(X)
            δθδY[i] = _δθδYᵢ(i, (1/β₁), Ω, V, U, sY, sX, σxy)
        end
        δβ₀δX = zeros(AbstractFloat, length(X))
        Threads.@threads for i in eachindex(X)
            δβ₀δX[i] = _δβ₀δXᵢ(i, Ȳ, (1/β₁), Ω, δθδX, δθδβ₁)
        end
        δβ₀δY = zeros(AbstractFloat, length(X))
        Threads.@threads for i in eachindex(X)
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
    return β₀, β₀SE, β₁, β₁SE, X_intercept, X_interceptSE, χ²ᵣ, pval, nX
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
    σxy::AbstractArray)
        sum(Ω.^2 .* (2β₁ .* (U .* V .* sX.^2 .- U.^2 .* σxy) .+ (U.^2 .* sY.^2 .- V.^2 .* sX.^2))) +
        4 * sum(Ω.^3 .* (σxy .- β₁ .* sX.^2) .* (
                β₁^2 .* (U .* V .* sX.^2 .- U.^2 .* σxy) .+
                β₁ .* (U.^2 .* sY.^2 .- V.^2 .* sX.^2) .-
                (U .* V .* sY.^2 .- V.^2 .* σxy))) +
        2 *
        sum(Ω.^2 .* (
            -(β₁^2) .* U .* sX.^2 .+
            2β₁ .* V .* sX.^2 .+
            U .* sY.^2 .- 2 .* V .* σxy)) *
        sum(Ω.^2 .* V .* (σxy .- β₁ .* sX.^2)) /
        sum(Ω) +
        2 *
        sum(Ω.^2 .* (
            -(β₁^2) .* V .* sX.^2 .+
            2(β₁^2) .* U .* σxy .-
            2β₁ .* U .* sY.^2 .+ V .* sY.^2)) *
        sum(Ω.^2 .* U .* (σxy .- β₁ .* sX.^2)) /
        sum(Ω)
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
    @simd for j in eachindex(Ω)
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
    @simd for j in eachindex(Ω)
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
    δθδβ₁::AbstractFloat
)
    return -β₁ * Ω[i] / sum(Ω) - X̄ * δθδX[i] / δθδβ₁
end

function _δβ₀δYᵢ(
    i::Integer,
    X̄::AbstractFloat,
    Ω::AbstractArray,
    δθδY::AbstractArray,
    δθδβ₁::AbstractFloat
)
    return Ω[i] / sum(Ω) - X̄ * δθδY[i] / δθδβ₁
end

# δβ₀X = (-β₁ .* Ω) ./ sum(Ω .- X̄ .* (δθX / δθδβ₁))
# δβ₀Y = Ω ./ sum(Ω .- X̄ .* (δθY / δθδβ₁))
