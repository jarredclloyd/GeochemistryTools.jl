#= Preamble
This file contains functions for errors-in-variables linear regression.
=#

export yorkfit

#Caller functions

"""
    yorkfit(df::AbstractDataFrame; [se_level_in::Int=2, se_level_out::Int=2, se_type::AbstractString="abs",
        initial::Any=nothing])

Compute line of best fit using York errors-in-variables linear regression algorithm.

Input df as a DataFrame of 4 of 5 columns wide with column order (X, σx, y, σy, [ρxy]).

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

# References
York D. et al. 2004 "Unified equations for the slope, intercept, and standard errors of the best straight line", *American
Journal of Physics*, 72(3), doi:https://doi.org/10.1119/1.1632486.

"""
function yorkfit(
    df::AbstractDataFrame;
    se_level_in::Int = 2,
    # se_level_out::Int = 2,
    se_type::AbstractString = "abs",
    initial::Any = nothing,
)
    dfCols::Int = ncol(df)
    dfRows::Int = nrow(df)
    if initial !== nothing
        if isa(initial, String) == true && haskey(dict_sr87_sr86i, initial) == true
            initial = deepcopy(dict_sr87_sr86i[initial])
        elseif length(initial) == 1 && isa(initial, AbstractFloat)
            initial = [1e-15, 1e-15, initial, 0.01]
        elseif length(initial) == 2 && isa(initial[1], AbstractFloat) && isa(initial[2], AbstractFloat)
            initial = [1e-15, 1e-15, initial[1], initial[2]]
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
            yfit = _eivlr_york(
                df[!, 1], df[!, 2] ./ se_level_in, df[!, 3], df[!, 4] ./ se_level_in, df[!, 5]
            )
        elseif dfCols == 4
            yfit = _eivlr_york(
                df[!, 1], df[!, 2] ./ se_level_in, df[!, 3], df[!, 4] ./ se_level_in
            )
        else
            throw(ArgumentError("Column width is not equal to 4 or 5. Some data is missing."))
        end
    elseif occursin('r', lowercase.(se_type)) == true ||
        occursin("rel", lowercase.(se_type)) == true ||
        occursin("relative", lowercase.(se_type)) == true
        if dfCols == 5
            yfit = _eivlr_york(
                df[!, 1],
                (df[!, 2] .* df[!, 1]) ./ se_level_in,
                df[!, 3],
                (df[!, 4] .* df[!, 3]) ./ se_level_in,
                df[!, 5],
            )
        elseif dfCols == 4
            yfit = _eivlr_york(
                df[!, 1], (df[!, 2] .* df[!, 1]) ./ se_level_in, df[!, 3], (df[!, 4] .* df[!, 3]) ./ se_level_in
            )
        else
            throw(ArgumentError("Column width is not equal to 4 or 5. Some data is missing."))
        end
    end
    if initial !== nothing
        popat!(df, dfRows + 1)
    end
    return yfit
end

function _eivlr_york(x::Vector{<:Real}, σx::Vector{<:Real}, y::Vector{<:Real}, σy::Vector{<:Real}, ρxy::Union{Nothing, Vector{<:Real}} = nothing)
    if _check_equal_length(x, σx, y, σy) != true
        throw(ArgumentError("The length of x, σx, y and σy must be the same."))
    end
    𝑁::Int = length(x)
    if ρxy === nothing
        ρxy::Vector{Float64} = zeros(𝑁)
    elseif length(ρxy) !== 𝑁
        throw(ArgumentError("The length of ρxy must be the same as x, σx, y and σy."))
    end

    # initial slope and intercept from OLS
    β₀::Float64, β₁::Float64 = hcat(ones(𝑁), x) \ y
    βₑ::Float64 = β₁

    # weights
    ωx::Vector{Float64} = @. 1.0 / σx ^ 2
    ωy::Vector{Float64} = @. 1.0 / σy ^ 2

    # initial fit via York method
    α::Vector{Float64} = sqrt.(ωx .* ωy)
    Ω::Vector{Float64} = @. ωx * ωy / (ωx + β₁^2 * ωy - 2 * β₁ * ρxy * α)
    x̄::Float64 = sum(Ω .* x) / sum(Ω)
    ȳ::Float64 = sum(Ω .* y) / sum(Ω)
    u::Vector{Float64} = x .- x̄
    v::Vector{Float64} = y .- ȳ
    βᵢ::Vector{Float64} = @. Ω * (u / ωy + β₁ * v / ωx - (β₁ * u + v) * ρxy / α)
    β₁ = sum(Ω .* βᵢ .* v) / sum(Ω .* βᵢ .* u)

    n_iterations::Int = 1
    # iterative solve via York method
    while (βₑ / β₁ - 1)^2 > eps() && n_iterations < 1e6
        Ω = @. ωx * ωy / (ωx + β₁^2 * ωy - 2 * β₁ * ρxy * α)
        x̄ = sum(Ω .* x) / sum(Ω)
        ȳ = sum(Ω .* y) / sum(Ω)
        u = x .- x̄
        v = y .- ȳ
        βᵢ = @. Ω * (u / ωy + β₁ * v / ωx - (β₁ * u + v) * ρxy / α)
        βₑ = β₁
        β₁ = sum(Ω .* βᵢ .* v) / sum(Ω .* βᵢ .* u)
        n_iterations += 1
    end
    β₀ = ȳ - β₁ * x̄
    xᵢ::Vector{Float64} = x̄ .+ βᵢ
    x̄ = sum(Ω .* xᵢ) / sum(Ω)
    u = xᵢ .- x̄
    β₁SE::Float64 = √(1 / sum(Ω .* u .^ 2))
    β₀SE::Float64 = √(1 / sum(Ω) + (x̄ * β₁SE)^2)
    σᵦ₁ᵦ₀::Float64 = - x̄ * β₁SE^2
    χ²::Float64 = sum(Ω .* (y .- β₁ .* x .- β₀) .^ 2)
    ν::Int = 𝑁 > 2 ? 𝑁 - 2 : 1
    χ²ᵣ::Float64 = χ² / ν
    pval::Float64 = ccdf(Chisq(ν), χ²)
    x_intercept::Float64 = -β₀ / β₁
    x_intercept_se::Float64 = sqrt((β₀SE / β₀)^2 + (β₁SE / β₁)^2 - 2 * σᵦ₁ᵦ₀ / (β₀ * β₁))
    return York(β₀, β₀SE, β₁, β₁SE, x_intercept, x_intercept_se, χ²ᵣ, pval, σᵦ₁ᵦ₀, 𝑁, x̄, ȳ)
end
