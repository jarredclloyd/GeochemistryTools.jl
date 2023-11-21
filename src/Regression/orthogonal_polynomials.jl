#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-09-08
Edited: 2023-09-19

This source file contains functions to compute orthogonal polynomial fits (up to pₙ(5)) and
their uncertainties. These are based on the equations provided in Bevington & Robinson 2003,
and Anenburg & Williams (2022).

Bevington, PR & Robinson, DK (2003) 'Data reduction and error analysis for the physical
sciences', 3rd ed., McGraw-Hill, ISBN: 9780072472271

Anenburg, M & Williams, MJ (2022) 'Quantifying the Tetrad Effect, Shape Components, and
Ce–Eu–Gd Anomalies in Rare Earth Element Patterns', *Mathematical Geosciences*, 54(1):47–70.
https://doi.org/10.1007/s11004-021-09959-5

=#
# function exports
export fit_orthogonal
export poly_orthogonal, poly_confidenceband, poly_predictionband, poly_standarderror

# stucts and base extensions
struct OrthogonalPolynomial <: LinearRegression
    lambda::AbstractVector
    lambda_se::AbstractArray
    beta::Real
    gamma::AbstractVector
    delta::AbstractVector
    epsilon::AbstractVector
    variance_covariance::AbstractArray
    order::AbstractVector
    r_squared::AbstractVector
    rmse::AbstractVector
    chi_squared::AbstractVector
    reduced_chi_squared::AbstractVector
    akaike_information_criteria::AbstractVector
    bayesian_information_criteria::AbstractVector
    n_observations::Integer
end

function Base.show(io::IOContext, fit::OrthogonalPolynomial)
    println(io, "λ₀: $(round(fit.lambda[1], sigdigits = 5))")
    println(io, "λ₁: $(round(fit.lambda[2], sigdigits = 5))")
    println(io, "λ₂: $(round(fit.lambda[3], sigdigits = 5))")
    println(io, "λ₃: $(round(fit.lambda[4], sigdigits = 5))")
    println(io, "λ₄: $(round(fit.lambda[5], sigdigits = 5))")
end

# call functions

function fit_orthogonal(
    df::AbstractDataFrame,
    x_name::Symbol,
    y_name::Symbol;
    y_weights::Union{Nothing,Symbol} = nothing,
    weight_by::AbstractString = "abs",
    rm_outlier::Bool = false,
)
    if y_weights !== nothing
        return _orthogonal_LSQ(
            df[!, x_name],
            df[!, y_name];
            y_weights = df[!, y_weights],
            weight_by = weight_by,
            rm_outlier = rm_outlier,
        )
    else
        return _orthogonal_LSQ(
            df[!, x_name],
            df[!, y_name];
            weight_by = weight_by,
            rm_outlier = rm_outlier,
        )
    end
end

function fit_orthogonal(
    A::AbstractArray;
    errors::Bool = false,
    weight_by::AbstractString = "abs",
    rm_outlier::Bool = false,
)
    if errors === false
        return _orthogonal_LSQ(
            A[:, 1],
            A[:, 2];
            weight_by = weight_by,
            rm_outlier = rm_outlier,
        )
    elseif errors === true
        return _orthogonal_LSQ(
            A[:, 1],
            A[:, 2];
            y_weights = A[:, 3],
            weight_by = weight_by,
            rm_outlier = rm_outlier,
        )
    end
end

function poly_orthogonal(x::AbstractVector, fit::OrthogonalPolynomial, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    return _poly_orthogonal(
        x,
        fit.lambda,
        fit.beta,
        fit.gamma,
        fit.delta,
        fit.epsilon,
        order,
    )
end

function poly_standarderror(
    x,
    fit::OrthogonalPolynomial,
    order::Integer;
    se_level::Integer = 2,
)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    X = _design_matrix(x, fit, order)
    VarΛX = fit.variance_covariance[1:(order + 1), 1:(order + 1)]
    return vec(
        sqrt.((fit.rmse[order + 1]^2) .* sum(X .* (X * VarΛX); dims = 2)) .* se_level,
    )
end

function poly_confidenceband(
    x,
    fit::OrthogonalPolynomial,
    order::Integer;
    ci_level::AbstractFloat = 0.95,
)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    tvalue = cquantile(TDist(length(x) - order), (1 - ci_level) / 2)
    X = _design_matrix(x, fit, order)
    VarΛX = fit.variance_covariance[1:(order + 1), 1:(order + 1)]
    return vec(sqrt.((fit.rmse[order + 1]^2) .* sum(X .* (X * VarΛX); dims = 2)) .* tvalue)
end

function poly_predictionband(
    x,
    fit::OrthogonalPolynomial,
    order::Integer;
    ci_level::AbstractFloat = 0.95,
)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    tvalue = cquantile(TDist(length(x) - order), (1 - ci_level) / 2)
    X = _design_matrix(x, fit, order)
    VarΛX = fit.variance_covariance[1:(order + 1), 1:(order + 1)]
    return vec(
        sqrt.((fit.rmse[order + 1]^2) .* sum(1 .+ X .* (X * VarΛX); dims = 2)) .* tvalue,
    )
end

# primary calculation function
function _orthogonal_LSQ(
    x::AbstractVector,
    y::AbstractVector;
    y_weights::Union{Nothing,AbstractArray} = nothing,
    weight_by::AbstractString = "abs",
    rm_outlier::Bool = false,
)
    𝑁 = length(x)
    β = _beta_orthogonal(x)
    γ = _gamma_orthogonal(x)
    δ = _delta_orthogonal(x)
    ϵ = _epsion_orthogonal(x)
    order = [0, 1, 2, 3, 4]
    X = hcat(
        repeat([1.0], 𝑁),
        (x .- β),
        (x .- γ[1]) .* (x .- γ[2]),
        (x .- δ[1]) .* (x .- δ[2]) .* (x .- δ[3]),
        (x .- ϵ[1]) .* (x .- ϵ[2]) .* (x .- ϵ[3]) .* (x .- ϵ[4]),
    )
    if y_weights === nothing
        ω = repeat([1.0], length(y))
    elseif occursin("abs", lowercase(weight_by)) === true
        ω = y_weights
    elseif occursin("rel", lowercase(weight_by)) == true
        ω = y_weights ./ y
    else
        throw(
            ArgumentError(
                "Value of 'weight_by' is unrecognised. String should contain either 'rel' or 'abs'.",
            ),
        )
    end
    ω = ω ./ mean(ω)
    ω = 1 ./ ω .^ 2
    Ω = Diagonal(ω)
    if rm_outlier === true
        VarΛX = inv(transpose(X) * (Ω) * X)
        Λ = VarΛX * transpose(X) * Ω * y
        Ĥ = X * VarΛX * transpose(X) * (Ω)
        mse4 = (transpose((y .- (X * Λ))) * Ω * (y .- (X * Λ))) / (𝑁 .- 5)
        studentised_residuals = (y .- (X * Λ)) ./ (sqrt.(mse4 .* (1 .- diag(Ĥ))))
        X = @view X[Not(studentised_residuals .>= 3), :]
        y = @view y[Not(studentised_residuals .>= 3)]
        ω = @view ω[Not(studentised_residuals .>= 3)]
    end
    Ω = Diagonal(ω)
    VarΛX = inv(transpose(X) * (Ω) * X)
    Λ = VarΛX * transpose(X) * Ω * y
    rss = zeros(5)
    @simd for i ∈ eachindex(order)
        @inbounds rss[i] =
            transpose((y .- (X[:, 1:i] * Λ[1:i]))) * Ω * (y .- (X[:, 1:i] * Λ[1:i]))
    end
    mse = rss ./ (𝑁 .- (order .+ 1))
    Λ_SE = zeros(5, 5)
    @simd for i ∈ eachindex(order)
        @inbounds Λ_SE[1:i, i] = sqrt.(diag(VarΛX[1:i, 1:i] * (mse[i])))
    end
    tss = transpose((y .- mean(y))) * Ω * (y .- mean(y))
    rmse = sqrt.(mse)
    R² = 1 .- (rss ./ (tss))
    for i ∈ eachindex(R²)
        if R²[i] < 0
            R²[i] = 0
        else
            R²[i] = _olkin_pratt.(R²[i], 𝑁, order[i] + 1)
        end
    end
    AIC = zeros(5)
    BIC = zeros(5)
    for i ∈ eachindex(order)
        AIC[i] = _akaike_information_criteria(rss[i], 𝑁, order[i])
        BIC[i] = _bayesian_information_criteria(rss[i], 𝑁, order[i])
    end
    return OrthogonalPolynomial(
        Λ,
        Λ_SE,
        β,
        γ,
        δ,
        ϵ,
        VarΛX,
        order,
        R²,
        rmse,
        rss,
        mse,
        AIC,
        BIC,
        𝑁,
    )
end

# polynomial functions
function _poly_orthogonal(
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

# functions for parameter calculations
function _beta_orthogonal(x::AbstractVector)
    return 1 / length(x) * sum(x)
end

function _gamma_orthogonal(x::AbstractVector)
    vieta = [-sum(x) length(x); -(sum(x .^ 2)) sum(x)] \ [-sum(x .^ 2); -sum(x .^ 3)]
    return real(roots(Polynomial([vieta[2], -vieta[1], 1]); permute = false, scale = false))
end

function _delta_orthogonal(x::AbstractVector)
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

function _epsion_orthogonal(x::AbstractVector)
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

function _design_matrix(x::AbstractVector, fit::OrthogonalPolynomial, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    if order == 0
        X = repeat([1.0], length(x))
    elseif order == 1
        X = hcat(repeat([1.0], length(x)), (x .- fit.beta))
    elseif order == 2
        X = hcat(
            repeat([1.0], length(x)),
            (x .- fit.beta),
            (x .- fit.gamma[1]) .* (x .- fit.gamma[2]),
        )
    elseif order == 3
        X = hcat(
            repeat([1.0], length(x)),
            (x .- fit.beta),
            (x .- fit.gamma[1]) .* (x .- fit.gamma[2]),
            (x .- fit.delta[1]) .* (x .- fit.delta[2]) .* (x .- fit.delta[3]),
        )
    elseif order == 4
        X = hcat(
            repeat([1.0], length(x)),
            (x .- fit.beta),
            (x .- fit.gamma[1]) .* (x .- fit.gamma[2]),
            (x .- fit.delta[1]) .* (x .- fit.delta[2]) .* (x .- fit.delta[3]),
            (x .- fit.epsilon[1]) .* (x .- fit.epsilon[2]) .* (x .- fit.epsilon[3]) .*
            (x .- fit.epsilon[4]),
        )
    end
    return X
end
