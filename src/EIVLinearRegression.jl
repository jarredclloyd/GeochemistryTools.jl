#= Preamble
This file contains functions for error-in-variable linear regression.
=#

export york, yorkfit, deming, demingfit, jackknife

#Caller functions
"""
    demingfit(df::DataFrame)

    Compute line of best fit using Deming errors-in-variables regression algorithm.
"""
function demingfit(df::DataFrame)
    β₀::Float64, β₁::Float64 = deming(df[!, 1]::Vector{Float64}, df[!, 2]::Vector{Float64})
    if isnan(β₀) || isnan(β₁)
        println("Solution not computed!")
    else
        β₀SE::Float64, β₁SE::Float64, n = jackknife(df, "deming")
    end
    return β₀, β₀SE, β₁, β₁SE, n
end

"""
    yorkfit(df::DataFrame; [SElevel::Int=2, SEtype::String="abs", SrInitial::Any=nothing])

Compute line of best fit using York errors-in-variables regression algorithm.

Input df as a DataFrame of 4 of 5 columns wide with column order (X, sX, Y, sY, [ρXY]).

# Optional arguments 
SElevel: standard error level, input as an integer. \n
SEtype: standard error type as a string of value "abs" OR "rel". \n
SrInitial: ⁸⁷Sr/⁸⁶Srᵢ value. Can be input as a string key from `Dict_SrInitial`, as a single numeric value, or as a vector of the initial and its standard error (same SElevel as input data). E.g. SrInitial = "MDCInv", SrInitial = 0.72, OR SrInitial = [0.72, 0.01].

For a full list of available keys in `Dict_SrInitial`, type `keys(Dict_SrInitial)`

# Example
```julia-repl
julia> york(df, 2, "abs"; SrInitial = "MDCInv")

```

# Reference
York et al. 2004 doi:https://doi.org/10.1119/1.1632486.
"""
function yorkfit(df::DataFrame; SElevel::Int=2, SEtype::String="abs", SrInitial::Any=nothing)
    dfCols::Int64 = ncol(df)
    if SrInitial != nothing
        if isa(SrInitial, String) == true && haskey(Dict_SrInitial, SrInitial) == true
            push!(df, Dict_SrInitial[SrInitial])
        elseif length(SrInitial) == 1 && isa(SrInitial, Number)
            push!(df, [1e-10, 1e-18, SrInitial, 0.01])
        elseif length(SrInitial) == 2 && isa(SrInitial[1], Number) && isa(SrInitial[2], Number)
            push!(df, [1e-10, 1e-18, SrInitial[1], SrInitial[2]])
        elseif length(SrInitial) >= 3
            error("Only the initial Sr ratio value and its uncertainty is required.")
        else
            error("""
            SrInitial variable is malformed. 
            Use either a string that matches a fixed value or reference material name: 
                e.g. "MDC", "MDCinv", or "SolarSystem" 
            Or input a single numeric value or two-length numeric vector: 
                e.g. 0.72 or [0.7200, 0.0010] 
             
            For a full list of available reference material keys: keys(Dict_SrInitial)
            """)
        end
    end
    if SEtype == "abs"
        if dfCols == 5
            β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX = york(df[:, 1], df[:, 2] ./ SElevel, df[:, 3], df[:, 4] ./ SElevel,
                df[:, 5])
        elseif dfCols == 4
            β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX = york(df[:, 1], df[:, 2] ./ SElevel, df[:, 3], df[:, 4] ./ SElevel)
        else
            println("Column width is not equal to 4 or 5. Some data is missing.")
        end
    elseif SEtype == "rel"
        if dfCols == 5
            β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX = york(df[:, 1], (df[:, 2] .* df[:, 1]) ./ SElevel, df[:, 3],
                (df[:, 4] .* df[:, 3]) ./ SElevel, df[:, 5])
        elseif dfCols == 4
            β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX = york(df[:, 1], (df[:, 2] .* df[:, 1]) ./ SElevel, df[:, 3],
                (df[:, 4] .* df[:, 3]) ./ SElevel)
        else
            println("Column width is not equal to 4 or 5. Some data is missing.")
        end
    end
    return β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX
end

#Base functions
function deming(X, Y)
    nX::Int = length(X)
    if nX ≠ length(Y)
        println("Lengths of data vectors are not equal. Please ensure all data (X & Y) are of equal length.")
    else
        λ::Float64 = (var(Y)) / (var(X))
        X̄::Float64 = mean(X)
        Ȳ::Float64 = mean(Y)
        Sxx::Float64 = sum((X .- X̄) .^ 2)
        Syy::Float64 = sum((Y .- Ȳ) .^ 2)
        Sxy::Float64 = sum((X .- X̄) .* (Y .- Ȳ))
        β₁::Float64 = (Syy - λ * Sxx + √((Syy - λ * Sxx)^2 + 4 * λ * Sxy^2)) / (2 * Sxy)
        β₀::Float64 = Ȳ - β₁ * X̄
    end
    return β₀, β₁
end

#jackknife deming
function jackknife(df::DataFrame, method::String="deming")
    nX::Int = nrow(df)
    if method == "deming"
        jackknifed_estimates::Matrix = zeros(Float64, nX, 2)
        @simd for j_position ∈ 1:nX
            holder = popat!(df, j_position)
            β₀est::Float64, β₁est::Float64 = deming(df[!, 1], df[!, 2])
            jackknifed_estimates[j_position, 1] = β₀est
            jackknifed_estimates[j_position, 2] = β₁est
            insert!(df, j_position, holder)
        end
        β₀SE::Float64 = sem(jackknifed_estimates[:, 1])
        β₁SE::Float64 = sem(jackknifed_estimates[:, 2])
        return β₀SE, β₁SE, nX
    end
end

function york(X::Vector{Float64}, sX::Vector{Float64}, Y::Vector{Float64}, sY::Vector{Float64}, ρXY=[0.0, 0.0])
    nX::Int = length(X)
    if length(ρXY) ≠ nX
        ρXY::Vector{Float64} = zeros(nX)
    end
    β₀::Float64, β₁::Float64 = coeffs(Polynomials.fit(X, Y, 1))
    βₑ::Float64 = β₁
    ωXᵢ::Vector{Float64} = 1 ./ sX .^ 2
    ωYᵢ::Vector{Float64} = 1 ./ sY .^ 2
    α::Vector{Float64} = .√(ωXᵢ .* ωYᵢ)
    Ω::Vector{Float64} = ωXᵢ .* ωYᵢ ./ (ωXᵢ .+ β₁^2 .* ωYᵢ .- 2 .* β₁ .* ρXY .* α)
    X̄::Float64 = sum(Ω .* X) / sum(Ω)
    Ȳ::Float64 = sum(Ω .* Y) / sum(Ω)
    U::Vector{Float64} = X .- X̄
    V::Vector{Float64} = Y .- Ȳ
    βᵢ::Vector{Float64} = Ω .* (U ./ ωYᵢ + β₁ * V ./ ωXᵢ - (β₁ * U + V) .* ρXY ./ α)
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
    x̄::Float64 = sum(Ω .* X) / sum(Ω)
    υ::Vector{Float64} = X .- x̄
    β₁SE::Float64 = √(1 / sum(Ω .* υ .^ 2))
    β₀SE::Float64 = √(1 / sum(Ω) + (x̄ * β₁SE)^2)
    σᵦ₁ᵦ₀::Float64 = -x̄ * β₁SE^2
    χ²::Float64 = sum(Ω .* (Y .- β₁ .* X .- β₀) .^ 2)
    ν::Int = nX - 2
    χ²ᵣ::Float64 = χ² / ν
    pval::Float64 = ccdf(Chisq(ν), χ²)
    return β₀, β₀SE, β₁, β₁SE, χ²ᵣ, pval, σᵦ₁ᵦ₀, nX
end