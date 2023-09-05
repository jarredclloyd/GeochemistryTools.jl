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
