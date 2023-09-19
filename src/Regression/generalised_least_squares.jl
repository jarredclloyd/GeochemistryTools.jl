#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-09-16
Edited: 2023-09-16

This source file contains functions to perform generalised least squares.

β = (XᵀΩX)⁻¹XᵀΩy
Cov(β|X) | C = (XᵀΩX)⁻¹
H = X(XᵀC⁻¹X)⁻¹XᵀX⁻¹
(error sum of squares) ess = (y - Xβ)ᵀΩ(y - Xβ)
(total sum of squares) tss = (y - ȳ)ᵀΩ(y - ȳ)
MSE = ess / 𝑁
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
    y_weights::Union{Nothing,AbstractVector} = nothing,
    weight_by::AbstractString = "abs",
)
    order = abs(order)
    if length(x) != length(y)
        throw(ArgumentError("Lengths of vectors 'x' and 'y' are not equal!"))
    elseif y_weights !== nothing && length(y) != length(y_weights)
        throw(ArgumentError("Lengths of vectors 'y' and 'y_weights' are not equal!"))
    end
    if y_weights === nothing
        ω = repeat([1.0], length(y))
    elseif occursin("abs", lowercase(weight_by)) === true
        ω = y_weights
    elseif occursin("rel", lowercase(weight_by)) == true
        ω = y_weights ./ y
    else
        throw(
            ArgumentError(
                "Value of 'weight_by' is unrecognised. String should contain either 'rel'
                or 'abs'.",
            ),
        )
    end
    ω = ω ./ median(ω)
    ω = 1 ./ ω .^ 2
    Ω = Diagonal(ω)
    X = _design_matrix(x, order)
    C = inv(transpose(X) * Ω * X)
    β = inv(transpose(X) * Ω * X) * transpose(X) * Ω * y
    ess = transpose(y .- X * β) * Ω * (y .- X * β)
    β_SE = sqrt.(diag(abs.((C) * (ess / (length(x) - order)))))
    tss = transpose((y .- mean(y))) * Ω * (y .- mean(y))
    mse = ess / (length(x) - order)
    R² = 1 - (ess / tss)
    if R² < 0
        R² = 0
    end
    R² = _olkin_pratt(R², length(x), order)
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
    order = length(fit.β) - 1
    X = _design_matrix(x, order)
    return X * fit.β[1:(order + 1)]
end

function _polyCI(
    x::AbstractVector,
    fit::GeneralisedLeastSquares;
    ci_level::AbstractFloat = 0.95,
)
    order = length(fit.β) - 1
    tvalue = cquantile(TDist(fit.𝑁 - order), (1 - ci_level) / 2)
    X = _design_matrix(x, order)
    return vec(sqrt.((fit.rmse^2) .* sum(X .* (X * fit.VarβX); dims = 2))) .* tvalue
end

function _polyPI(
    x::AbstractVector,
    fit::GeneralisedLeastSquares;
    ci_level::AbstractFloat = 0.95,
)
    order = length(fit.β) - 1
    tvalue = cquantile(TDist(fit.𝑁 - order), (1 - ci_level) / 2)
    X = _design_matrix(x, order)
    return vec(
        sqrt.((fit.rmse^2) .* sum(1 .+ X .* (X * fit.VarβX) .* (fit.rmse^2); dims = 2)),
    ) .* tvalue
end
