function _olkin_pratt(R²::AbstractFloat, 𝑁::Integer, predictors::Integer)
    z = 1 - R²
    _₂F₁ = HypergeometricFunctions._₂F₁positive(1, 1, (𝑁 - predictors + 1) / 2, z)
    return 1 - ((𝑁 - 3) / (𝑁 - predictors + 1) * z) * _₂F₁
end

function _chi_squared_reduced(χ²::Real, 𝑁::Integer, predictors::Integer)
    return χ² / (𝑁 - predictors)
end


function _bayesian_information_criteria(χ²::AbstractFloat, 𝑁::Integer, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    return χ² + order * log(𝑁)
end

function _chi_squared(
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
        r = y .- _poly_orthogonal(x, λ, β, γ, δ, ϵ, order)
    elseif order == 2
        r = y .- _poly_orthogonal(x, λ, β, γ, δ, ϵ, order)
    elseif order == 3
        r = y .- _poly_orthogonal(x, λ, β, γ, δ, ϵ, order)
    elseif order == 4
        r = y .- _poly_orthogonal(x, λ, β, γ, δ, ϵ, order)
    end
    χ² = transpose(r) * Ω * r
    return χ²
end