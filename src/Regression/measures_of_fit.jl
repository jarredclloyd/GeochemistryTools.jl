function _olkin_pratt(R²::AbstractFloat, 𝑛::Integer, predictors::Integer)
    z = 1 - R²
    c = (𝑛 - predictors + 1) / 2
    _₂F₁value = HypergeometricFunctions._₂F₁positive(1, 1, c, z)
    return 1 - ((𝑛 - 3) / (𝑛 - predictors - 1)) * z * _₂F₁value
end

function _chi_squared_reduced(χ²::Real, 𝑛::Integer, predictors::Integer)
    return χ² / (𝑛 - predictors)
end


function _bayesian_information_criteria(rss::Real, 𝑛::Integer, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    k = order + 2
    return 𝑛 * log(rss / 𝑛) + k * log(𝑛) + 𝑛 * log(2π) + 𝑛
end


function _akaike_information_criteria(rss::Real, 𝑛::Integer, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    k = order + 2
    return 𝑛 * log(rss / 𝑛) + 2 * k + ((2 * k * (k + 1)) / (𝑛 - k - 1))
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

function _reduced_chi_squared_ci(dof::Integer, confidence_level::AbstractFloat = 0.95)
    lower_χ²ᵣ = cquantile(Chisq(dof), 1 - (1 - confidence_level) / 2) / dof
    upper_χ²ᵣ = cquantile(Chisq(dof), (1 - confidence_level) / 2) / dof
    return (lower_χ²ᵣ, upper_χ²ᵣ)
end
