#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Created: 2023-09-08
Edited: 2023-11-29

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
    lambda::Vector{Float64}
    lambda_se::SparseMatrixCSC{Float64}
    beta::Float64
    gamma::Vector{Float64}
    delta::Vector{Float64}
    epsilon::Vector{Float64}
    variance_covariance::Symmetric{Float64,Matrix{Float64}}
    order::Vector{Integer}
    r_squared::Vector{Float64}
    rmse::Vector{Float64}
    chi_squared::Vector{Float64}
    reduced_chi_squared::Vector{Float64}
    akaike_information_criteria::Vector{Float64}
    bayesian_information_criteria::Vector{Float64}
    n_observations::Integer
end

function Base.show(io::IOContext, fit::OrthogonalPolynomial)
    println(io, "λ₀: $(round(fit.lambda[1], sigdigits = 5))")
    println(io, "λ₁: $(round(fit.lambda[2], sigdigits = 5))")
    println(io, "λ₂: $(round(fit.lambda[3], sigdigits = 5))")
    println(io, "λ₃: $(round(fit.lambda[4], sigdigits = 5))")
    return println(io, "λ₄: $(round(fit.lambda[5], sigdigits = 5))")
end

# call functions
"""

    fit_orthogonal(df::AbstractDataFrame,
    x_name::Symbol,
    y_name::Symbol;
    [y_weights::Union{Nothing,Symbol} = nothing,
    weight_by::AbstractString = "abs",
    rm_outlier::Bool = false,
    verbose::Bool = false])

Compute an orthogonal polynomial that represent some X and Y data.

Input df as a DataFrame of 4 of 5 columns wide with column order (X, sX, Y, sY, [ρXY]).

# Keywords
- `y_weights::Union{Nothing,Symbol}`: Weights for y values (e.g. absolute uncertainties).
- 'weight_by::AbstractString': Weight pre-scaling, values of "rel" or "abs" (default) are
    accepted. If "rel" transforms weights to relative weights.
- 'rm_outlier::Bool': When set to true, will remove outliers (studentised residuals ≥ 3,
    based on fit with minimum akaike information criteria value).
- 'verbose::Bool': When set to true will print the number of outliers determined during N passes.

# References

Bevington, PR & Robinson, DK (2003) 'Data reduction and error analysis for the physical
sciences', 3rd ed., McGraw-Hill, ISBN: 9780072472271

Anenburg, M & Williams, MJ (2022) 'Quantifying the Tetrad Effect, Shape Components, and
Ce–Eu–Gd Anomalies in Rare Earth Element Patterns', *Mathematical Geosciences*, 54(1):47–70.
https://doi.org/10.1007/s11004-021-09959-5

Akaike, H (1974) 'A new look at the statistical model identification',
*IEEE Transactions on Automatic Control*, 19(6):716–723.
https://doi.org/10.1109/TAC.1974.1100705

Karch, J (2020) 'Improving on Adjusted R-Squared', *Collabra: Psychology*, 6(1):45.
https://doi.org/10.1525/collabra.343

"""
function fit_orthogonal(
    df::AbstractDataFrame,
    x_name::Symbol,
    y_name::Symbol;
    y_weights::Union{Nothing,Symbol} = nothing,
    weight_by::AbstractString = "abs",
    rm_outlier::Bool = false,
    verbose::Bool = false
)
    if y_weights !== nothing
        return _orthogonal_LSQ(
            df[!, x_name],
            df[!, y_name];
            y_weights = df[!, y_weights],
            weight_by = weight_by,
            rm_outlier = rm_outlier,
            verbose = verbose
        )
    else
        return _orthogonal_LSQ(
            df[!, x_name],
            df[!, y_name];
            weight_by = weight_by,
            rm_outlier = rm_outlier,
            verbose = verbose
        )
    end
end

function fit_orthogonal(
    A::AbstractArray;
    errors::Bool = false,
    weight_by::AbstractString = "abs",
    rm_outlier::Bool = false,
    verbose::Bool = false
)
    if errors === false
        return _orthogonal_LSQ(
            A[:, 1],
            A[:, 2];
            weight_by = weight_by,
            rm_outlier = rm_outlier,
            verbose = verbose
        )
    elseif errors === true
        return _orthogonal_LSQ(
            A[:, 1],
            A[:, 2];
            y_weights = A[:, 3],
            weight_by = weight_by,
            rm_outlier = rm_outlier,
            verbose = verbose
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
    X::Matrix{Float64} = _design_matrix(x, fit, order)
    VarΛX = view(fit.variance_covariance, 1:(order + 1), 1:(order + 1))
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
    X::Matrix{Float64} = _design_matrix(x, fit, order)
    VarΛX = view(fit.variance_covariance, 1:(order + 1), 1:(order + 1))
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
    X::Matrix{Float64} = _design_matrix(x, fit, order)
    VarΛX = view(fit.variance_covariance, 1:(order + 1), 1:(order + 1))
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
    verbose::Bool = false,
)
    𝑁::Integer = length(x)
    x_sums::Vector{Float64} = Vector{Float64}(undef, 7)
    @simd for i ∈ eachindex(x_sums)
        x_sums[i] = sum(x .^ i)
    end
    β::Float64 = _beta_orthogonal(𝑁, x_sums)
    γ::Vector{Float64} = _gamma_orthogonal(𝑁, x_sums)
    δ::Vector{Float64} = _delta_orthogonal(𝑁, x_sums)
    ϵ::Vector{Float64} = _epsilon_orthogonal(𝑁, x_sums)
    order::Vector{Int64} = [0, 1, 2, 3, 4]
    X::Matrix{Float64} = hcat(
        fill(1.0, 𝑁),
        (x .- β),
        (x .- γ[1]) .* (x .- γ[2]),
        (x .- δ[1]) .* (x .- δ[2]) .* (x .- δ[3]),
        (x .- ϵ[1]) .* (x .- ϵ[2]) .* (x .- ϵ[3]) .* (x .- ϵ[4]),
    )
     if y_weights === nothing
        ω::Vector{Float64} = fill(1.0, length(y))
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
    Ω::Diagonal{Float64,Vector{Float64}} = Diagonal(1 ./ (ω ./ mean(ω)) .^ 2)
    Xᵀ::Transpose{Float64,Matrix{Float64}} = transpose(X)
    rss::Vector{Float64} = Vector{Float64}(undef, 5)
    AIC::Vector{Float64} = Vector{Float64}(undef, 5)
    VarΛX::Symmetric{Float64,Matrix{Float64}} = Symmetric(inv(Xᵀ * (Ω) * X))
    Λ::Vector{Float64} = VarΛX * Xᵀ * Ω * y
    @inbounds @simd for i ∈ eachindex(order)
        residuals::Vector{Float64} = (y .- (view(X, :, 1:i) * Λ[1:i]))
        rss[i] =
            transpose(residuals) * Ω * (residuals)
    end
    AIC = _akaike_information_criteria.(rss, 𝑁, order)
    if rm_outlier === true
        𝑁prev::Integer = 0
        n_iterations::Integer = 0
        n_outliers::Integer = 0
        while 𝑁prev - 𝑁 != 0 && n_iterations ≤ 10
            n_iterations += 1
            minAIC::Int64 = findmin(AIC)[2]
            Xvar::Matrix{Float64} = view(VarΛX, 1:minAIC, 1:minAIC) * view(Xᵀ, 1:minAIC, :)
            leverage::Vector{Float64} = Vector{Float64}(undef, size(X, 1))
            Threads.@threads for i ∈ axes(X, 1)
                @inbounds leverage[i] = sum(view(X, i, 1:minAIC) .* view(Xvar, :, i))
            end
            leverage .*= ω
            studentised_residuals::Vector{Float64} =
                y .- (view(X, :, 1:minAIC) * Λ[1:minAIC]) # 3 allocs
            mse::Vector{Float64} = rss ./ (𝑁 .- (order .+ 1))
            studentised_residuals ./= @.(sqrt(mse[minAIC] * (1 - leverage)))
            outlier_inds::Vector{Int64} = findall(studentised_residuals .≥ 3)
            n_outliers += length(outlier_inds)
            if n_outliers > 0
                X = view(X, Not(outlier_inds), :) # high allocs
                y = y[Not(outlier_inds)] # high allocs
                ω = ω[Not(outlier_inds)] # high allocs
                Xᵀ = transpose(X)
                Ω = Diagonal(1 ./ (ω ./ mean(ω)) .^ 2)

                VarΛX = Symmetric(inv(Xᵀ * (Ω) * X))
                Λ = VarΛX * Xᵀ * Ω * y
                @inbounds @simd for i ∈ eachindex(order)
                    residuals = (y .- (view(X, :, 1:i) * Λ[1:i]))
                    rss[i] = transpose(residuals) * Ω * (residuals)
                end
                AIC = _akaike_information_criteria.(rss, 𝑁, order)
            end
            𝑁prev = 𝑁
            𝑁 = size(X, 1)
        end
        if verbose == true
            println(
                "Determined $n_outliers $(n_outliers == 1 ?  "outlier" : "outliers") for current fit in $n_iterations $(n_iterations == 1 ?  "pass" : "passes")",
            )
        end
    end
    mse = rss ./ (𝑁 .- (order .+ 1))
    Λ_SE::AbstractMatrix{Float64} = zeros(Float64, 5, 5)
    @inbounds for i ∈ eachindex(order)
        Λ_SE[1:i, i] = sqrt.(diag(view(VarΛX, 1:i, 1:i) * (mse[i])))
    end
    sparse(Λ_SE)
    tss::Float64 = transpose((y .- mean(y))) * Ω * (y .- mean(y))
    rmse::Vector{Float64} = sqrt.(mse)
    R²::Vector{Float64} = 1 .- (rss ./ (tss))
    @inbounds for i ∈ eachindex(R²)
        if R²[i] < 0
            R²[i] = 0
        else
            R²[i] = _olkin_pratt(R²[i], 𝑁, order[i] + 1)
        end
    end
    BIC::Vector{Float64} = Vector{Float64}(undef, 5)
    BIC = _bayesian_information_criteria.(rss, 𝑁, order)
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
function _beta_orthogonal(N::Integer, sums::AbstractVector)
    return 1 / N * sums[1]
end

function _gamma_orthogonal(N::Integer, sums::AbstractVector)
    vieta::Vector{Float64} = [-sums[1] N; -sums[2] sums[1]] \ [-sums[2]; -sums[3]]
    return real(PolynomialRoots.roots(([vieta[2], -vieta[1], 1])))
end

function _delta_orthogonal(N::Integer, sums::AbstractVector)
    vieta::Vector{Float64} =
        [
            -sums[2] sums[1] -N
            -sums[3] sums[2] -sums[1]
            -sums[4] sums[3] -sums[2]
        ] \ [-sums[3]; -sums[4]; -sums[5]]
    return real(PolynomialRoots.roots(([-vieta[3], vieta[2], -vieta[1], 1])))
end
function _epsilon_orthogonal(N::Integer, sums::AbstractVector)
    vieta::Vector{Float64} =
        [
            -sums[3] sums[2] -sums[1] N
            -sums[4] sums[3] -sums[2] sums[1]
            -sums[5] sums[4] -sums[3] sums[2]
            -sums[6] sums[5] -sums[4] sums[3]
        ] \ [-sums[4]; -sums[5]; -sums[6]; -sums[7]]
    return real(PolynomialRoots.roots(([vieta[4], -vieta[3], vieta[2], -vieta[1], 1])))
end

function _design_matrix(x::AbstractVector, fit::OrthogonalPolynomial, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    if order == 0
        X::Matrix{Float64} = repeat([1.0], length(x))
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

function _squaredmahalanobis(n, hii)
    return (n - 1) * (hii - 1 / n)
end
