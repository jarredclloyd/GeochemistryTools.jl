#= Preamble
This file contains functions for errors-in-variables linear regression.
=#

export york, yorkfit, deming, demingfit, jackknife

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
    while (βₑ / β₁ - 1)^2 > 1e-15 && n_iterations < 1000
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
    x̄ = sum(Ω .* X) / sum(Ω)
    υ::AbstractArray{AbstractFloat} = X .- x̄
    β₁SE::AbstractFloat = √(1 / sum(Ω .* υ .^ 2))
    β₀SE::AbstractFloat = √(1 / sum(Ω) + (x̄ * β₁SE)^2)
    σᵦ₁ᵦ₀::AbstractFloat = - x̄ * β₁SE^2
    χ²::AbstractFloat = sum(Ω .* (Y .- β₁ .* X .- β₀) .^ 2)
    ν::Int = nX - 2
    χ²ᵣ::AbstractFloat = χ² / ν
    pval::AbstractFloat = ccdf(Chisq(ν), χ²)
    return β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX
end

#=
function MonteCarloEIV()

end
=#
