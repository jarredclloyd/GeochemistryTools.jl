#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-08-22
Edited: 2023-08-22

This source file contains functions to compute a line of best fit and its standard errors using the 'Deming'
errors-in-variables regression algorithm of Mahon (1996) corrected by Stephan and Trappitsch (2023).

=#
# function exports
export eivlr_fit
#Caller functions
"""
    demingfit(df::AbstractDataFrame)

Compute line of best fit using Deming errors-in-variables regression algorithm.

Input df as a DataFrame of 2 columns wide with order (x, y). Standard errors for β₀ and β₁.
"""
function demingfit(df::AbstractDataFrame)
    β₀::Float64, β₁::Float64 = deming(df[!, 1], df[!, 2])
    if isnan(β₀) || isnan(β₁)
        println("Solution not computed!")
    else
        β₀SE::Float64, β₁SE::Float64, n = jackknife(df, "deming")
    end
    return β₀, β₀SE, β₁, β₁SE, n
end

function deming(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    𝑁::Int = length(x)
    if 𝑁 !== length(y)
        throw(
            ArgumentError("Lengths of data vectors are not equal. Please ensure all data (x & y) are of equal length.")
        )
    else
        λ::Float64 = (var(y)) / (var(x))
        x̄::Float64 = mean(x)
        ȳ::Float64 = mean(y)
        Sxx::Float64 = sum((x .- x̄) .^ 2)
        Syy::Float64 = sum((y .- ȳ) .^ 2)
        Sxy::Float64 = sum((x .- x̄) .* (y .- ȳ))
        β₁::Float64 = (Syy - λ * Sxx + √((Syy - λ * Sxx)^2 + 4 * λ * Sxy^2)) / (2 * Sxy)
        β₀::Float64 = ȳ - β₁ * x̄
    end
    return β₀, β₁
end

#jackknife deming
function jackknife(df::AbstractDataFrame, method::AbstractString="deming")
    𝑁::Int = nrow(df)
    if method == "deming"
        jackknifed_estimates::AbstractVector{Float64} = zeros(Float64, 𝑁, 2)
        @threads for j_position in 1:𝑁
            β₀est::Float64, β₁est::Float64 = deming(df[Not(j_position), 1], df[Not(j_position), 2])
            jackknifed_estimates[j_position, 1] = β₀est
            jackknifed_estimates[j_position, 2] = β₁est
        end
        β₀SE::Float64 = sem(jackknifed_estimates[:, 1])
        β₁SE::Float64 = sem(jackknifed_estimates[:, 2])
        return β₀SE, β₁SE, 𝑁
    end
end
