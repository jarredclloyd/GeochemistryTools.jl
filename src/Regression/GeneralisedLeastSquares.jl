#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-09-16
Edited: 2023-09-16

This source file contains functions to perform generalised least squares.

β = (XᵀΩ⁻¹X)⁻¹XᵀΩ⁻¹y
Cov(β|X) | C = (XᵀΩ⁻¹X)⁻¹
H = X(XᵀC⁻¹X)⁻¹XᵀX⁻¹
(residual sum of squares) rss = (y - Xβ)ᵀΩ(y - Xβ)
(total sum of squares) tss = (y - ȳ)ᵀΩ(y - ȳ)
MSE = rss / (𝑁-k+1)
Model CI = yᵢ ± √(mse * row sum(X * (XC)) * t(1-α/2, n - p) where p is number of predictors (order)
Model PI = yᵢ ± √(mse .* row sum(1 + x * (XC)) * t(1-α/2, n - p)
(biased) R² = ess / tss
(Olkin-Pratt) R² = 1-((N - 3) / (N - p - 1)) * (1 - R²) * ₂F₁(1,1; (N - p +1 ) / 2; 1 - R²)

=#
# function exports

function _GLS(
    x::AbstractVector,
    y::AbstractVector,
    order::Integer;
    y_weights::Union{Nothing,AbstractVector}=nothing,
    weight_type::AbstractString="abs",
    adjust_r::Bool=false
)
    order = abs(order)
    if length(x) != length(y)
        throw(ArgumentError("Lengths of vectors 'x' and 'y' are not equal!"))
    elseif y_weights !== nothing && length(y) != length(y_weights)
        throw(ArgumentError("Lengths of vectors 'y' and 'y_weights' are not equal!"))
    end
    if y_weights === nothing
        ω = repeat([1.0], length(y))
    elseif occursin("abs", lowercase(weight_type)) === true
        ω = y_weights
    elseif occursin("rel", lowercase(weight_type)) == true
        ω = y_weights ./ y
    else
        throw(
            ArgumentError(
                "Value of 'weight_by' is unrecognised. String should contain either 'rel'
                or 'abs'.",
            ),
        )
    end
    ω = ω .^ 2
    Ω = Diagonal(ω)
    X = _design_matrix(x, order)
    C = inv(transpose(X) * inv(Ω) * X)
    β = inv(transpose(X) * inv(Ω) * X) * transpose(X) * inv(Ω) * y
    rss = transpose(y .- X * β) * inv(Ω) * (y .- X * β)
    β_SE = sqrt.(diag(abs.((C) * (rss / (length(x) - order)))))
    tss = transpose((y .- mean(y))) * inv(Ω) * (y .- mean(y))
    mse = rss / (length(x) - (order + 1))
    R² = 1 - (rss / tss)
    if R² ≤ Base.rtoldefault(Float64)
        R² = 0
    end
    if adjust_r == true
        R² = _olkin_pratt(R², length(x), order)
    end
    rmse = sqrt(mse)
    return GeneralisedLeastSquares(β, β_SE, C, R², rmse, length(x))
end

function _design_matrix(x, order)
    order = abs.(order)
    𝑁 = length(x)
    X = repeat([1.0], 𝑁)
    if order == 0
        return X
    else
        for i ∈ 1:order
            X = hcat(X, x .^ [i])
        end
    end
    return X
end

function _polyGLS(x::AbstractVector, fit::GeneralisedLeastSquares)
    order = length(fit.beta) - 1
    X = _design_matrix(x, order)
    return X * fit.beta[1:(order+1)]
end

function _polyCI(
    x::AbstractVector,
    fit::GeneralisedLeastSquares;
    ci_level::AbstractFloat=0.95,
)
    order = length(fit.beta) - 1
    tvalue = cquantile(TDist(fit.n_observations - order), (1 - ci_level) / 2)
    X = _design_matrix(x, order)
    return vec(sqrt.((fit.rmse^2) .* sum(X .* (X * fit.variance_covariance); dims=2))) .* tvalue
end

function _polyPI(
    x::AbstractVector,
    fit::GeneralisedLeastSquares;
    ci_level::AbstractFloat=0.95,
)
    order = length(fit.beta) - 1
    tvalue = cquantile(TDist(fit.n_observations - order), (1 - ci_level) / 2)
    X = _design_matrix(x, order)
    return vec(
        sqrt.((fit.rmse^2) .* sum(1 .+ X .* (X * fit.variance_covariance) .* (fit.rmse^2); dims=2)),
    ) .* tvalue
end
