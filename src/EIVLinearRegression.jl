#= Preamble
This file contains functions for error-in-variable linear regression.
=#

export york, yorkfit, deming, demingfit, jackknife, york2

#Caller functions
function demingfit(A::DataFrame)
    β₀::Float64, β₁::Float64 = deming(A[:, 1], A[:, 2])
    if isnan(β₀) || isnan(β₁)
        println("Solution not computed!")
    else
        β₀SE::Float64, β₁SE::Float64 = jackknife(A[:, 1], A[:, 2], "deming")
    end
    return β₀, β₀SE, β₁, β₁SE
end

function yorkfit(A::DataFrame)
    widthA::Int = ncol(A)
    if widthA == 5
        β₀, β₀SE, β₁, β₁SE, σᵦ₁ᵦ₀, χ², χ²ᵣ, pval = york(A[:, 1], A[:, 2], A[:, 3], A[:, 4], A[:, 5])
    elseif widthA == 4
        β₀, β₀SE, β₁, β₁SE, σᵦ₁ᵦ₀, χ², χ²ᵣ, pval = york(A[:, 1], A[:, 2], A[:, 3], A[:, 4])
    else
        println("Column width is not equal to 4 or 5. Some data is missing.")
    end
    return β₀, β₀SE, β₁, β₁SE, σᵦ₁ᵦ₀, χ², χ²ᵣ, pval
end

#Base functions
"""
    deming(X, Y)

    Compute line of best fit using Deming errors-in-variables regression algorithm.
"""
function deming(X::Vector{Float64}, Y::Vector{Float64})
    nX::Int = length(X)
    if nX ≠ length(Y)
        println("Lengths of data vectors are not equal. Please ensure all data (X & Y) are of equal length.")
    else
        β₀::Float64, β₁::Float64 = GLM.coef(lm(@formula(Y ~ X), DataFrame(X=X, Y=Y)))
        if isnan(β₀) || isnan(β₁)
            println("Cannot fit a straight line through these data.")
        else
            λ::Float64 = (var(Y)) / (var(X))
            X̄::Float64 = mean(X)
            Ȳ::Float64 = mean(Y)
            Sxx::Float64 = sum((X .- X̄) .^ 2)
            Syy::Float64 = sum((Y .- Ȳ) .^ 2)
            Sxy::Float64 = sum((X .- X̄) .* (Y .- Ȳ))
            β₁ = (Syy - λ * Sxx + √((Syy - λ * Sxx)^2 + 4 * λ * Sxy^2)) / (2 * Sxy)
            β₀ = Ȳ - β₁ * X̄
        end
        return β₀, β₁
    end
end

#jackknife deming
function jackknife(X::Vector{Float64}, Y::Vector{Float64}, method::String="deming")
    nX::Int = length(X)
    if method == "deming"
        jackknifed_estimates = zeros(Float64, nX, 2)
        @simd for j_position ∈ 1:nX
            jackknifedata::DataFrame = DataFrame(X=X, Y=Y)
            deleteat!(jackknifedata, j_position)
            β₀est::Float64, β₁est::Float64 = deming(jackknifedata[:, 1], jackknifedata[:, 2])
            jackknifed_estimates[j_position, :1] = β₀est
            jackknifed_estimates[j_position, :2] = β₁est
        end
        YintSE::Float64 = sem(jackknifed_estimates[:, 1])
        slopeSE::Float64 = sem(jackknifed_estimates[:, 2])
        return YintSE, slopeSE
    end
end

"""
    york(X, sX, Y, sY, ρXY)

Compute line of best fit using York errors-in-variables regression algorithm.

Based on York et al. 2004 doi:https://doi.org/10.1119/1.1632486

"""
function york(X::Vector{Float64}, sX::Vector{Float64}, Y::Vector{Float64}, sY::Vector{Float64}, ρXY::Vector{Float64}=[0.0, 0.0])
    nX::Int = length(X)
    if length(ρXY) ≠ nX
        ρXY::Vector{Float64} = repeat([ρₑ(sX, sY)], nX)
    end
    if nX ≠ length(Y) || nX ≠ length(sX) || nX ≠ length(sY) || nX ≠ length(ρXY)
        println("Lengths of data vectors are not equal. Please ensure all data (X,Y,sX,sY,ρXY) are of equal length.")
    else
        β₀::Float64, β₁::Float64 = GLM.coef(lm(@formula(Y ~ X), DataFrame(X=X, Y=Y)))
        if isnan(β₀) || isnan(β₁)
            println("Cannot fit a straight line through these data.")
        else
            βₑ::Float64 = β₁
            ωXᵢ::Vector{Float64} = 1 ./ sX .^2
            ωYᵢ::Vector{Float64} = 1 ./ sY .^2
            α::Vector{Float64} = .√(ωXᵢ .* ωYᵢ)
            Ω::Vector{Float64} = ωXᵢ .* ωYᵢ ./ (ωXᵢ .+ β₁^2 .* ωYᵢ .- 2 .* β₁ .* ρXY .* α)
            X̄::Float64 = sum(Ω .* X) / sum(Ω)
            Ȳ::Float64 = sum(Ω .* Y) / sum(Ω)
            U::Vector{Float64} = X .- X̄
            V::Vector{Float64} = Y .- Ȳ
            βᵢ::Vector{Float64} = Ω .* (U ./ ωYᵢ + β₁ * V ./ ωXᵢ - (β₁ * U + V) .* ρXY ./ α)
            β₁ = sum(Ω .* βᵢ .* V) / sum(Ω .* βᵢ .* U)
            n_iterations::Int = 1
            while (βₑ / β₁ - 1)^2 > 1e-24 && n_iterations < 1000
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
        end
        β₀ = Ȳ - β₁ * X̄
        x̄::Float64 = sum(Ω .* X) / sum(Ω)
        υ::Vector{Float64} = X .- x̄
        β₁SE::Float64 = √(1 / sum(Ω .* υ .^ 2))
        β₀SE::Float64 = √(1 / sum(Ω) + (x̄ * β₁SE)^2)
        σᵦ₁ᵦ₀::Float64 = -x̄ * β₁SE ^2
        χ²::Float64 = sum(Ω .*(Y .- β₁ .* X .- β₀) .^2)
        df::Int = nX - 2
        χ²ᵣ::Float64 = χ² / df
        pval::Float64 = ccdf(Chisq(df),χ²)
        return β₀, β₀SE, β₁, β₁SE, σᵦ₁ᵦ₀, χ², χ²ᵣ, pval
    end
end

#Helper functions

function ρₑ(sX, sY)
    ρₛ::Float64 = corspearman(sX, sY)
    ρₖ::Float64 = corkendall(sX, sY)
    ρₐ = [ρₛ, ρₖ]
    ρₑ::Float64 = mean(ρₐ)
end

