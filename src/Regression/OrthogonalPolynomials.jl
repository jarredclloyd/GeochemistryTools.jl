#= Preamble

Author: Jarred C Lloyd: https://github.com/jarredclloyd
Last edited: 2025-04-25

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

# structs and base extensions
struct OrthogonalPolynomial <: LinearRegression
    lambda::Union{Vector{AbstractFloat},Nothing}
    lambda_se::Union{AbstractMatrix,Nothing}
    beta::Union{AbstractFloat,Nothing}
    gamma::Union{Vector{AbstractFloat},Nothing}
    delta::Union{Vector{AbstractFloat},Nothing}
    epsilon::Union{Vector{AbstractFloat},Nothing}
    variance_covariance::Union{AbstractMatrix,Nothing}
    order::Union{Vector{Integer},Nothing}
    r_squared::Union{Vector{AbstractFloat},Nothing}
    OP_r_squared::Union{Vector{AbstractFloat},Nothing}
    rmse::Union{Vector{AbstractFloat},Nothing}
    nrmse::Union{Vector{AbstractFloat},Nothing}
    chi_squared::Union{Vector{AbstractFloat},Nothing}
    reduced_chi_squared::Union{Vector{AbstractFloat},Nothing}
    akaike_information_criteria::Union{Vector{AbstractFloat},Nothing}
    akaike_weights::Union{Vector{AbstractFloat},Nothing}
    bayesian_information_criteria::Union{Vector{AbstractFloat},Nothing}
    bayesian_weights::Union{Vector{AbstractFloat},Nothing}
    n_observations::Union{Integer,Nothing}
end

function Base.show(io::IOContext, fit::OrthogonalPolynomial)
    if isnothing(fit.n_observations)
        println(io, "no model was fitted to this data")
    else
        println(io, repeat("-", 80))
        println(io, "Orthogonal Polynomial Model fitted $(fit.n_observations) observations from input data")
        pretty_table(
            hcat(
                ["λ₀", "λ₁", "λ₂", "λ₃", "λ₄"],
                round.(Float64.(fit.lambda); sigdigits = 4),
                round.(fit.r_squared; sigdigits = 4),
                round.(fit.OP_r_squared; sigdigits = 4),
                round.(fit.reduced_chi_squared; sigdigits = 4),
                round.(fit.nrmse; sigdigits = 4),
                round.(fit.akaike_weights; sigdigits = 4),
                round.(fit.bayesian_weights; sigdigits = 4),
            );
            header = ["Lambda", "Value", "ρ²", "ρ²ₒₚ", "χ²ᵣ", "RMSE (normalised)", "AICc Weight", "BICc Weight"],
        )
        println(io, "\n Model Uncertainties (%1SE)")
        pretty_table(hcat(["λ₀", "λ₁", "λ₂", "λ₃", "λ₄"],
            round.(abs.(Float64.((fit.lambda_se) ./ fit.lambda) .* 100); sigdigits = 4));
            header = ["Lambda", "Φ₀(x)", "Φ₁(x)", "Φ₂(x)", "Φ₃(x)", "Φ₄(x)"],
        )
        println(io, repeat("-", 80))
    end
end

# call functions
"""
    fit_orthogonal(df::AbstractDataFrame,
    x_name::Symbol,
    y_name::Symbol;
    [y_weights::Union{Nothing,Symbol} = nothing,
    weight_type::AbstractString = "abs",
    rm_outlier::Bool = false,
    verbose::Bool = false])

Compute an orthogonal polynomial that represent some X and Y data.

Input df as a DataFrame and specify the x and y data column names as symbols.

# Keywords

  - `y_weights::Union{Nothing,Symbol}`: Weights for y values (e.g. absolute uncertainties).
  - `weight_type::AbstractString = "abs"`: Weight pre-scaling, values of "rel" or "abs" are
    accepted. If "abs" transforms weights to relative weights.
  - `rm_outlier::Bool = false`: When set to true, will remove outliers (studentised residuals ≥ 3,
    based on fit with minimum akaike information criteria value).
  - `verbose::Bool = false`: When set to true will print the number of outliers determined during N passes.

# Description
This function computes the orthogonal polynomial fit up to order five for some linear x-y data.
It can account for errors/uncertainties in y, and perform automated outlier removal.
If outlier removal is enabled, the algorithm will compute the standardised residuals for
the modelled polynomial with the minimised corrected akaike information criteria and remove
points with a standardised residual ≥ 3 until no outliers remain, or ten iterations have
been performed.

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

Burnham, KP & Anderson, DR (2002) 'Model selection and multimodel inference: A practical
information-theoretic approach', 2nd ed., Springer, ISBN: 978-0-387-95364-9
"""
function fit_orthogonal(
    df::AbstractDataFrame,
    x_name::Symbol,
    y_name::Symbol;
    y_weights::Union{Nothing,Symbol}=nothing,
    weight_type::AbstractString="abs",
    rm_outlier::Bool=false,
    verbose::Bool=false,
)
    if y_weights !== nothing
        return _orthogonal_LSQ(
            df[!, x_name],
            df[!, y_name];
            y_weights=df[!, y_weights],
            weight_type=weight_type,
            rm_outlier=rm_outlier,
            verbose=verbose,
        )
    else
        return _orthogonal_LSQ(
            df[!, x_name],
            df[!, y_name];
            weight_type=weight_type,
            rm_outlier=rm_outlier,
            verbose=verbose,
        )
    end
end

"""

    fit_orthogonal(A::AbstractArray;
    errors::Bool = false,
    weight_type::AbstractString = "abs",
    rm_outlier::Bool = false,
    verbose::Bool = false])

Compute an orthogonal polynomial that represent some X and Y data.

Input A as an Array of 4 of 5 columns wide with column order (X, sX, Y, sY, [ρXY]). If the

# Keywords

  - `y_weights::Union{Nothing,Symbol}`: Weights for y values (e.g. absolute uncertainties).
  - `weight_type::AbstractString = "abs"`: Weight pre-scaling, values of "rel" or "abs" are
    accepted. If "abs" transforms weights to relative weights.
  - `rm_outlier::Bool = false`: When set to true, will remove outliers (studentised residuals ≥ 3,
    based on fit with minimum akaike information criteria value).
  - `verbose::Bool = false`: When set to true will print the number of outliers determined during N passes.
  - `st_residual_tol::Real = 3.0`: Tolerance for studentised residuals for outlier detection. Y-values with a studentised residual ≥ this value will be considered an outlier.

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

Burnham, KP & Anderson, DR (2002) 'Model selection and multimodel inference: A practical
information-theoretic approach', 2nd ed., Springer, ISBN: 978-0-387-95364-9
"""
function fit_orthogonal(
    A::AbstractArray;
    errors::Bool=false,
    weight_type::AbstractString="rel",
    rm_outlier::Bool=false,
    verbose::Bool=false,
    st_residual_tol::Real=3.0
)
    if errors === false
        return _orthogonal_LSQ(
            A[:, 1],
            A[:, 2];
            weight_type=weight_type,
            rm_outlier=rm_outlier,
            verbose=verbose,
        )
    elseif errors === true
        return _orthogonal_LSQ(
            A[:, 1],
            A[:, 2];
            y_weights=A[:, 3],
            weight_type=weight_type,
            rm_outlier=rm_outlier,
            verbose=verbose,
        )
    end
end

function poly_orthogonal(x::Union{AbstractVector, Real}, fit::OrthogonalPolynomial, order::Integer)
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
    se_level::Integer=2,
)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    X = _design_matrix(x, fit, order)
    VarΛX = view(fit.variance_covariance, 1:(order+1), 1:(order+1))
    return vec(
        sqrt.(fit.reduced_chi_squared[order+1] .* sum(X .* (X * VarΛX); dims=2)) .* se_level,
    )
end

function poly_confidenceband(
    x,
    fit::OrthogonalPolynomial,
    order::Integer;
    ci_level::AbstractFloat=0.95,
)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    tvalue = cquantile(TDist(length(x) - order + 1), (1 - ci_level) / 2)
    X = _design_matrix(x, fit, order)
    VarΛX = view(fit.variance_covariance, 1:(order+1), 1:(order+1))
    return vec(
        sqrt.(fit.reduced_chi_squared[order+1] .* sum(X .* (X * VarΛX); dims=2)) .* tvalue
    )
end

function poly_predictionband(
    x,
    fit::OrthogonalPolynomial,
    order::Integer;
    ci_level::AbstractFloat=0.95,
)
    if order < 0
        throw(ArgumentError("Polynomial order must be positive"))
    end
    tvalue = cquantile(TDist(length(x) - order + 1), (1 - ci_level) / 2)
    X = _design_matrix(x, fit, order)
    VarΛX = view(fit.variance_covariance, 1:(order+1), 1:(order+1))
    return vec(
        sqrt.(fit.reduced_chi_squared[order+1] .* (1 .+ sum(X .* (X * VarΛX); dims=2))) .* tvalue,
    )
end

# primary calculation function
function _orthogonal_LSQ(
    x::AbstractVector,
    y::AbstractVector;
    y_weights::Union{Nothing,AbstractArray}=nothing,
    weight_type::AbstractString="abs",
    rm_outlier::Bool=false,
    verbose::Bool=false,
    st_residual_tol::Real=3.0
)
    finite_indices::Vector{Int64} = intersect(findall(isfinite, x), findall(isfinite, y))
    x::Vector{Float64x4} = Float64x4.(x[finite_indices])
    y::Vector{Float64x4} = Float64x4.(y[finite_indices])
    𝑁::Integer = length(x)
    if 𝑁 == length(y) && 𝑁 > 2
        x_sums::MVector{7, Float64x4} = MVector{7, Float64x4}(undef)
        @simd for i ∈ eachindex(x_sums)
            x_sums[i] = sum(x .^ i)
        end
        β::Float64x4 = _beta_orthogonal(𝑁, x_sums)
        γ::SVector{2,Float64x4} = _gamma_orthogonal(𝑁, x_sums)
        δ::SVector{3,Float64x4} = _delta_orthogonal(𝑁, x_sums)
        ϵ::SVector{4,Float64x4} = _epsilon_orthogonal(𝑁, x_sums)
        order::SVector{5,Int64} = [0, 1, 2, 3, 4]

        # Construct design matrix, minimises allocations
        X::Matrix{Float64x4} = Matrix{Float64x4}(undef, (𝑁, 5))
        X[:, 1] .= 1.0
        X[:, 2] = x .- β
        X[:, 3] = (x .- γ[1]) .* (x .- γ[2])
        X[:, 4] = (x .- δ[1]) .* (x .- δ[2]) .* (x .- δ[3])
        X[:, 5] = (x .- ϵ[1]) .* (x .- ϵ[2]) .* (x .- ϵ[3]) .* (x .- ϵ[4])
        if y_weights === nothing
            ω::Vector{Float64x4} = fill(1.0, length(y))
        elseif occursin("rel", lowercase(weight_type)) === true
            ω = y_weights[finite_indices] .* y[finite_indices]
        elseif occursin("abs", lowercase(weight_type)) == true
            ω = y_weights[finite_indices]
        else
            throw(
                ArgumentError(
                    "Value of 'weight_type' is unrecognised. String should contain either 'rel' or 'abs'.",
                ),
            )
        end
        Ω::Diagonal{Float64x4,Vector{Float64x4}} = Diagonal(ω .^ 2)
        X̃::Matrix{Float64x4} = exp(-0.5log(Ω)) * X
        ỹ::Vector{Float64x4} = exp(-0.5log(Ω)) * y

        if cond(X̃) ≤ 1e7
            if verbose
                println("Fitting via QR factorisation")
            end
            F = qr(X̃)
            Λ::Vector{Float64x4} = F \ ỹ
            VarΛX = Symmetric(inv(F.R) * transpose(inv(F.R)))
        else
            if verbose
                println("Fitting via SVD method")
            end
            F = svd(X̃)
            Λ = F.V * inv(Diagonal(F.S)) * transpose(F.U) * ỹ
            VarΛX = F.V * inv(Diagonal(F.S .^ 2)) * F.Vt
        end
        rss::MVector{5, Float64x4} = MVector{5, Float64x4}(undef)
        @simd for i ∈ eachindex(order)
            residuals::Vector{Float64x4} = (y .- (view(X, :, 1:i) * Λ[1:i]))
            rss[i] = transpose(residuals) * (residuals)
        end
        AIC::SVector{5,Float64x4} = _akaike_information_criteria.(rss, 𝑁, order)
        if rm_outlier === true
            𝑁prev::Int64 = 0
            n_iterations::Int64 = 0
            t_outliers::Int64 = 0
            while 𝑁prev - 𝑁 != 0 && n_iterations ≤ 10
                n_outliers::Int64 = 0
                n_iterations += 1
                minAIC::Int64 = findmin(AIC)[2]
                Xvar::Matrix{Float64x4} =
                    view(VarΛX, 1:minAIC, 1:minAIC) * view(X', 1:minAIC, :) * inv(Ω)
                leverage::Vector{Float64x4} =
                    Vector{Float64x4}(undef, size(X, 1))
                Threads.@threads for i ∈ axes(X, 1)
                    leverage[i] = sum(view(X, i, 1:minAIC) .* view(Xvar, :, i))
                end
                studentised_residuals::Vector{Float64x4} =
                    y .- (view(X, :, 1:minAIC) * Λ[1:minAIC])
                mse::MVector{5, Float64x4} = rss ./ (𝑁 .- (order .+ 1))
                studentised_residuals ./= @.(sqrt(mse[minAIC] * (1 - leverage)))
                outlier_inds::Vector{Int64} = findall(≥(st_residual_tol), abs.(studentised_residuals))
                n_outliers += length(outlier_inds)
                t_outliers += n_outliers
                if verbose
                    println("Pass number $n_iterations: Found $n_outliers outliers")
                end
                if n_outliers > 0
                    X = view(X, Not(outlier_inds), :) # high allocs
                    y = y[Not(outlier_inds)] # high allocs
                    ω = ω[Not(outlier_inds)] # high allocs
                    Ω = Diagonal(ω .^ 2)
                    X̃ = exp(-0.5log(Ω)) * X
                    ỹ = exp(-0.5log(Ω)) * y
                    if cond(X) ≤ 1e7
                        if verbose
                            println("Fitting via QR factorisation")
                        end
                        F = qr(X̃)
                        Λ = F \ ỹ
                        VarΛX = Symmetric(inv(F.R) * transpose(inv(F.R)))
                    else
                        if verbose
                            println("Fitting via SVD method")
                        end
                        F = svd(X̃)
                        Λ = F.V * inv(Diagonal(F.S)) * transpose(F.U) * ỹ
                        VarΛX = F.V * inv(Diagonal(F.S .^ 2)) .* F.Vt
                    end
                    @simd for i ∈ eachindex(order)
                        residuals = (y .- (view(X, :, 1:i) * Λ[1:i]))
                        rss[i] = transpose(residuals) * (residuals)
                    end
                    AIC = _akaike_information_criteria.(rss, 𝑁, order)
                end
                𝑁prev = 𝑁
                𝑁 = size(X, 1)
            end
            if verbose == true
                println(
                    "Determined $t_outliers $(t_outliers == 1 ?  "outlier" : "outliers") for current fit in $n_iterations $(n_iterations == 1 ?  "pass" : "passes")",
                )
            end
        end
        for i in eachindex(Λ)
            Λ[i] = abs(Λ[i]) ≤ Base.rtoldefault(Float64x4) ? 0.0 : Λ[i]
        end
        mse = rss ./ (𝑁 .- (order .+ 1))
        Λ_SE::AbstractMatrix{Float64x4} = zeros(Float64x4, 5, 5)
        for i ∈ eachindex(order)
            Λ_SE[1:i, i] = sqrt.(diag(view(VarΛX, 1:i, 1:i) * (mse[i])))
        end
        χ²::MVector{5, Float64x4} = MVector{5, Float64x4}(undef)
        @simd for i ∈ eachindex(order)
                residuals = (y .- (view(X, :, 1:i) * Λ[1:i]))
                χ²[i] = transpose(residuals) * inv(Ω) * (residuals)
        end
        χ²ᵣ::SVector{5,Float64} = _chi_squared_reduced.(χ², 𝑁, order  .+ 1)
        tss::Float64x4 = transpose((y .- mean(y))) * (y .- mean(y))
        rmse::SVector{5, Float64} = sqrt.(mse)
        nrmse::SVector{5, Float64} = rmse ./ (maximum(y) - minimum(y))
        R²::SVector{5, Float64} = 1 .- (rss ./ (tss))
        R²ₒₚ::SVector{5, Float64} = _olkin_pratt.(R², 𝑁, order .+ 1)
        AIC = _akaike_information_criteria.(rss, 𝑁, order)
        BIC::SVector{5,Float64} = _bayesian_information_criteria.(rss, 𝑁, order)
        BICw::SVector{5, Float64} =
            exp.(-0.5 .* (BIC .- minimum(BIC))) ./ sum(exp.(-0.5 .* (BIC .- minimum(BIC))))
        AICw::SVector{5, Float64} =
            exp.(-0.5 .* (AIC .- minimum(AIC))) ./ sum(exp.(-0.5 .* (AIC .- minimum(AIC))))
        return OrthogonalPolynomial(
            Λ,
            UpperTriangular(Λ_SE),
            β,
            γ,
            δ,
            ϵ,
            VarΛX,
            order,
            R²,
            R²ₒₚ,
            rmse,
            nrmse,
            χ²,
            χ²ᵣ,
            AIC,
            AICw,
            BIC,
            BICw,
            𝑁,
        )
    else
        throw(error(("Unable to fit data as there are less than three values")))
        return OrthogonalPolynomial(
            fill(nothing, length(fieldnames(OrthogonalPolynomial)))...,
        )
    end
end

# polynomial functions
function _poly_orthogonal(
    x::Union{AbstractVector, Real},
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
        @. (λ[1] + 0 * x)
    elseif order == 1
        @. (λ[1] + λ[2] * (x - β))
    elseif order == 2
        @. (λ[1] + λ[2] * (x - β) + λ[3] * ((x - γ[1]) * (x - γ[2])))
    elseif order == 3
        @. (
            λ[1] +
            λ[2] * (x - β) +
            λ[3] * ((x - γ[1]) * (x - γ[2])) +
            λ[4] * ((x - δ[1]) * (x - δ[2]) * (x - δ[3]))
        )
    elseif order == 4
        @. (
            λ[1] +
            λ[2] * (x - β) +
            λ[3] * ((x - γ[1]) * (x - γ[2])) +
            λ[4] * ((x - δ[1]) * (x - δ[2]) * (x - δ[3])) +
            λ[5] * ((x - ϵ[1]) * (x - ϵ[2]) * (x - ϵ[3]) * (x - ϵ[4]))
        )
    end
end

# functions for parameter calculations
function _beta_orthogonal(N::Integer, sums::AbstractVector)
    return Float64x4(1) / Float64x4(N) * sums[1]
end

function _gamma_orthogonal(N::Integer, sums::AbstractVector)
    vieta::Vector{Float64x4} =
        qr([-sums[1] N; -sums[2] sums[1]]) \ [-sums[2]; -sums[3]]
    return real(PolynomialRoots.roots([vieta[2], -vieta[1], 1]; polish=true, epsilon=eps(Float64x4)))
end

function _delta_orthogonal(N::Integer, sums::AbstractVector)
    vieta::Vector{Float64x4} =
        qr([
            -sums[2] sums[1] -N
            -sums[3] sums[2] -sums[1]
            -sums[4] sums[3] -sums[2]
        ]) \ [-sums[3]; -sums[4]; -sums[5]]
    return real(PolynomialRoots.roots([-vieta[3], vieta[2], -vieta[1], 1]; polish=true, epsilon=eps(Float64x4)))
end
function _epsilon_orthogonal(N::Integer, sums::AbstractVector)
    vieta::Vector{Float64x4} =
        qr([
            -sums[3] sums[2] -sums[1] N
            -sums[4] sums[3] -sums[2] sums[1]
            -sums[5] sums[4] -sums[3] sums[2]
            -sums[6] sums[5] -sums[4] sums[3]
        ]) \ [-sums[4]; -sums[5]; -sums[6]; -sums[7]]
    return real(PolynomialRoots.roots([vieta[4], -vieta[3], vieta[2], -vieta[1], 1]; polish=true, epsilon=eps(Float64x4)))
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

function _squaredmahalanobis(n, hii)
    return (n - 1) * (hii - 1 / n)
end
