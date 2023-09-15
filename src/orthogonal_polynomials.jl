#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-09-08
Edited: 2023-09-09

This source file contains functions to compute orthogonal polynomial fits (up to pₙ(5)) and their uncertainties. These
are based on the equations provided in Bevington & Robinson 2003, and Anenburg & Williams (2022).

Bevington, PR & Robinson, DK (2003) 'Data reduction and error analysis for the physical sciences', 3rd ed,
McGraw-Hill, ISBN: 9780072472271
Anenburg, M & Williams, MJ (2022) 'Quantifying the Tetrad Effect, Shape Components, and Ce–Eu–Gd Anomalies in Rare
Earth Element Patterns', *Mathematical Geosciences*, 54(1):47–70. https://doi.org/10.1007/s11004-021-09959-5

=#
# function exports
export fit_orth
export polyλ, polyCI

# stucts and base extensions
struct OrthogonalPoly
    λ::AbstractVector
    λ_SE::AbstractVector
    β::Real
    γ::AbstractVector
    δ::AbstractVector
    ϵ::AbstractVector
    Var::AbstractArray
    order::AbstractVector
    χ²::AbstractVector
    χ²ᵣ::AbstractVector
    BIC::AbstractVector
end

function Base.show(io::IOContext, fit::OrthogonalPoly)
    println(
        io,
        "λ₀: $(round(fit.λ[1], sigdigits = 5)) ± $(round(fit.λ_SE[1], sigdigits = 5)); p(1): χ²ᵣ = $(round(fit.χ²ᵣ[1], sigdigits = 5)), BIC = $(round(fit.BIC[1], sigdigits = 5))",
    )
    println(
        io,
        "λ₁: $(round(fit.λ[2], sigdigits = 5)) ± $(round(fit.λ_SE[2], sigdigits = 5)); p(2): χ²ᵣ = $(round(fit.χ²ᵣ[2], sigdigits = 5)), BIC = $(round(fit.BIC[2], sigdigits = 5))",
    )
    println(
        io,
        "λ₂: $(round(fit.λ[3], sigdigits = 5)) ± $(round(fit.λ_SE[3], sigdigits = 5)); p(3): χ²ᵣ = $(round(fit.χ²ᵣ[3], sigdigits = 5)), BIC = $(round(fit.BIC[3], sigdigits = 5))",
    )
    println(
        io,
        "λ₃: $(round(fit.λ[4], sigdigits = 5)) ± $(round(fit.λ_SE[4], sigdigits = 5)); p(4): χ²ᵣ = $(round(fit.χ²ᵣ[4], sigdigits = 5)), BIC = $(round(fit.BIC[4], sigdigits = 5))",
    )
    return println(
        io,
        "λ₄: $(round(fit.λ[5], sigdigits = 5)) ± $(round(fit.λ_SE[5], sigdigits = 5)); p(5): χ²ᵣ = $(round(fit.χ²ᵣ[5], sigdigits = 5)), BIC = $(round(fit.BIC[5], sigdigits = 5))",
    )
end

# call functions
function fit_orth(
    df::AbstractDataFrame,
    x_name::Symbol,
    y_name::Symbol,
    yσ_name::Union{Nothing,Symbol} = nothing;
    weight::AbstractString = "abs",
    weight_transform::AbstractString = "inv",
)
    if yσ_name !== nothing
        return _orthogonal_LSQ(
            df[!, x_name],
            df[!, y_name],
            df[!, yσ_name];
            weight = weight,
            weight_transform = weight_transform,
        )
    else
        return _orthogonal_LSQ(
            df[!, x_name],
            df[!, y_name];
            weight = weight,
            weight_transform = weight_transform,
        )
    end
end

function fit_orth(
    A::AbstractArray;
    errors::Bool = false,
    weight::AbstractString = "abs",
    weight_transform::AbstractString = "inv",
)
    if errors === false
        return _orthogonal_LSQ(
            A[:, 1],
            A[:, 2];
            weight = weight,
            weight_transform = weight_transform,
        )
    elseif errors === true
        return _orthogonal_LSQ(
            A[:, 1],
            A[:, 2],
            A[:, 3];
            weight = weight,
            weight_transform = weight_transform,
        )
    end
end

function polyλ(x::AbstractVector, fλ::OrthogonalPoly, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    return _polyλ(x, fλ.λ, fλ.β, fλ.γ, fλ.δ, fλ.ϵ, order)
end

function polyCI(x, fλ::OrthogonalPoly, order::Integer; CIlevel::AbstractFloat = 0.95)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    tvalue = cquantile(TDist(length(x) - order), (1 - CIlevel) / 2)
    X = _design_matrix(x, fλ, order)
    Var = fλ.Var[1:(order + 1), 1:(order + 1)]
    return vec(sqrt.(sum(X .* (X * Var); dims = 2)) .* tvalue)
end

# primary calculation function
function _orthogonal_LSQ(
    x::AbstractVector,
    y::AbstractVector,
    yσ::Union{Nothing,AbstractArray} = nothing;
    weight::AbstractString = "abs",
    weight_transform::AbstractString = "inv",
)
    𝑁 = length(x)
    β = _β(x)
    γ = _γ(x)
    δ = _δ(x)
    ϵ = _ϵ(x)
    order = [0, 1, 2, 3, 4]
    X = hcat(
        repeat([1.0], length(x)),
        (x .- β),
        (x .- γ[1]) .* (x .- γ[2]),
        (x .- δ[1]) .* (x .- δ[2]) .* (x .- δ[3]),
        (x .- ϵ[1]) .* (x .- ϵ[2]) .* (x .- ϵ[3]) .* (x .- ϵ[4]),
    )
    if yσ !== nothing && occursin("rel", lowercase(weight)) == true
        ω = (yσ ./ y)
    elseif yσ !== nothing && occursin("abs", lowercase(weight)) == true
        ω = yσ
    else
        ω = repeat([1], length(x))
    end
    if weight_transform == "log"
        ω = log.(ω)
    else
        ω = ω .^ 2
    end
    Ω = inv(Diagonal(ω))
    Λ = inv(transpose(X) * Ω * X) * transpose(X) * Ω * y
    Var = inv(transpose(X) * Ω * X)
    ŷ = X * Λ
    r = y .- ŷ
    # varΒ = diag((inv(transpose(X) * Ω * X)) * ((transpose(r) * Ω * r)/(N - order[i]))) # the diagonal of (sum of squares / df) * VarCov Matrix)
    # println(varΒ)
    Λ_SE = sqrt.([Var[1, 1], Var[2, 2], Var[3, 3], Var[4, 4], Var[5, 5]])
    χ² = zeros(5)
    χ²ᵣ = zeros(5)
    BIC = zeros(5)
    for i ∈ eachindex(order)
        χ²[i] = _χ²(x, y, Ω, Λ, β, γ, δ, ϵ, order[i])
        χ²ᵣ[i] = _χ²ᵣ(𝑁, χ²[i], order[i])
        BIC[i] = _bayesian_information_criteria(χ²[i], 𝑁, order[i])
    end
    return OrthogonalPoly(Λ, Λ_SE, β, γ, δ, ϵ, Var, order, χ², χ²ᵣ, BIC)
end

# polynomial functions
function _polyλ(
    x::AbstractVector,
    λ::AbstractVector,
    β::AbstractFloat,
    γ::AbstractVector,
    δ::AbstractVector,
    ϵ::AbstractVector,
    order::Integer,
)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    if order == 0
        λ[1] .+ 0 .* x
    elseif order == 1
        λ[1] .+ λ[2] .* (x .- β)
    elseif order == 2
        λ[1] .+ λ[2] .* (x .- β) .+ λ[3] .* ((x .- γ[1]) .* (x .- γ[2]))
    elseif order == 3
        λ[1] .+ λ[2] .* (x .- β) .+ λ[3] .* ((x .- γ[1]) .* (x .- γ[2])) .+
        λ[4] .* ((x .- δ[1]) .* (x .- δ[2]) .* (x .- δ[3]))
    elseif order == 4
        λ[1] .+ λ[2] .* (x .- β) .+ λ[3] .* ((x .- γ[1]) .* (x .- γ[2])) .+
        λ[4] .* ((x .- δ[1]) .* (x .- δ[2]) .* (x .- δ[3])) .+
        λ[5] .* ((x .- ϵ[1]) .* (x .- ϵ[2]) .* (x .- ϵ[3]) .* (x .- ϵ[4]))
    end
end

# statistical functions
function _bayesian_information_criteria(χ²::AbstractFloat, 𝑁::Integer, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    return χ² + order * log(𝑁)
end

function _χ²(
    x::AbstractVector,
    y::AbstractVector,
    Ω::AbstractArray,
    λ::AbstractVector,
    β::AbstractFloat,
    γ::AbstractVector,
    δ::AbstractVector,
    ϵ::AbstractVector,
    order::Integer,
)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    if order == 0
        r = y .- λ[1]
    elseif order == 1
        r = y .- _polyλ(x, λ, β, γ, δ, ϵ, order)
    elseif order == 2
        r = y .- _polyλ(x, λ, β, γ, δ, ϵ, order)
    elseif order == 3
        r = y .- _polyλ(x, λ, β, γ, δ, ϵ, order)
    elseif order == 4
        r = y .- _polyλ(x, λ, β, γ, δ, ϵ, order)
    end
    χ² = transpose(r) * Ω * r
    return χ²
end

function _χ²ᵣ(𝑁, χ², order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    else
        χ²ᵣ = χ² / (𝑁 - order)
    end
    return χ²ᵣ
end

# functions for parameter calculations
function _β(x)
    return 1 / length(x) * sum(x)
end

function _γ(x)
    vieta = [-sum(x) length(x); -(sum(x .^ 2)) sum(x)] \ [-sum(x .^ 2); -sum(x .^ 3)]
    return real(roots(Polynomial([vieta[2], -vieta[1], 1]); permute = false, scale = false))
end

function _δ(x)
    vieta =
        [
            -sum(x .^ 2) sum(x) -length(x)
            -sum(x .^ 3) sum(x .^ 2) -sum(x)
            -sum(x .^ 4) sum(x .^ 3) -sum(x .^ 2)
        ] \ [-sum(x .^ 3); -sum(x .^ 4); -sum(x .^ 5)]
    return real(
        roots(
            Polynomial([-vieta[3], vieta[2], -vieta[1], 1]);
            permute = false,
            scale = false,
        ),
    )
end

function _ϵ(x)
    vieta =
        [
            -sum(x .^ 3) sum(x .^ 2) -sum(x) length(x)
            -sum(x .^ 4) sum(x .^ 3) -sum(x .^ 2) sum(x)
            -sum(x .^ 5) sum(x .^ 4) -sum(x .^ 3) sum(x .^ 2)
            -sum(x .^ 6) sum(x .^ 5) -sum(x .^ 4) sum(x .^ 3)
        ] \ [-sum(x .^ 4); -sum(x .^ 5); -sum(x .^ 6); -sum(x .^ 7)]
    return real(
        roots(
            Polynomial([vieta[4], -vieta[3], vieta[2], -vieta[1], 1]);
            permute = false,
            scale = false,
        ),
    )
end

function _design_matrix(x, fλ, order)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    if order == 0
        X = repeat([1.0], length(x))
    elseif order == 1
        X = hcat(repeat([1.0], length(x)), (x .- fλ.β))
    elseif order == 2
        X = hcat(repeat([1.0], length(x)), (x .- fλ.β), (x .- fλ.γ[1]) .* (x .- fλ.γ[2]))
    elseif order == 3
        X = hcat(
            repeat([1.0], length(x)),
            (x .- fλ.β),
            (x .- fλ.γ[1]) .* (x .- fλ.γ[2]),
            (x .- fλ.δ[1]) .* (x .- fλ.δ[2]) .* (x .- fλ.δ[3]),
        )
    elseif order == 4
        X = hcat(
            repeat([1.0], length(x)),
            (x .- fλ.β),
            (x .- fλ.γ[1]) .* (x .- fλ.γ[2]),
            (x .- fλ.δ[1]) .* (x .- fλ.δ[2]) .* (x .- fλ.δ[3]),
            (x .- fλ.ϵ[1]) .* (x .- fλ.ϵ[2]) .* (x .- fλ.ϵ[3]) .* (x .- fλ.ϵ[4]),
        )
    end
    return X
end
