export affine_fit, affine_ci, affine_pi

abstract type LinearRegression <: Any end
abstract type ErrorsInVariablesRegression <: LinearRegression end

Base.getindex(LR::LinearRegression, i::Integer) = getfield(LR, i)
Base.keys(LR::LinearRegression) = Base.keys(LR)

struct York <: ErrorsInVariablesRegression
    beta0::Real
    beta0_se::Real
    beta1::Real
    beta1_se::Real
    x_intercept::Real
    x_intercept_se::Real
    reduced_chi_squared::Real
    p_value::Real
    covariance_beta::Real
    n_observations::Integer
    x_bar::Real
    y_bar::Real
end

function Base.show(io::IOContext, F::York)
    println(io, "βₒ: $(round(F.beta0, sigdigits = 5)) ± $(round(F.beta0_se, sigdigits = 5))")
    println(io, "β₁: $(round(F.beta1, sigdigits = 5)) ± $(round(F.beta1_se, sigdigits = 5))")
    println(io, "x_intercept: $(round(F.x_intercept, sigdigits = 5)) ± $(round(F.x_intercept_se, sigdigits = 5))")
    println(io, "χ²ᵣ (MSWD): $(round(F.reduced_chi_squared, sigdigits = 3)); p-value: $(round(F.p_value, sigdigits = 5))")
    println(io, "no. observations: $(F.n_observations)")
end

struct MahonNonFixed <: ErrorsInVariablesRegression
    beta0::Real
    beta0_se::Real
    beta1::Real
    beta1_se::Real
    x_intercept::Real
    x_intercept_se::Real
    reduced_chi_squared::Real
    p_value::Real
    n_observations::Integer
end

function Base.show(io::IOContext, F::MahonNonFixed)
    println(io, "βₒ: $(round(F.beta0, sigdigits = 5)) ± $(round(F.beta0_se, sigdigits = 5))")
    println(io, "β₁: $(round(F.beta1, sigdigits = 5)) ± $(round(F.beta1_se, sigdigits = 5))")
    println(io, "x_intercept: $(round(F.x_intercept, sigdigits = 5)) ± $(round(F.x_intercept_se, sigdigits = 5))")
    println(io, "χ²ᵣ (MSWD): $(round(F.reduced_chi_squared, sigdigits = 3)); p-value: $(round(F.p_value, sigdigits = 5))")
    println(io, "no. observations: $(F.n_observations)")
end

struct MahonFixed <: ErrorsInVariablesRegression
    beta0::Real
    beta0_se::Real
    beta1::Real
    beta1_se::Real
    x_intercept::Real
    x_intercept_se::Real
    reduced_chi_squared::Real
    p_value::Real
    n_observations::Integer
    fixed_point::Tuple{Real,Real,Real,Real}
end

function Base.show(io::IOContext, F::MahonFixed)
    println(
        io,
        "βₒ: $(round(F.beta0, sigdigits = 5)) ± $(round(F.beta0_se, sigdigits = 5))",
    )
    println(
        io,
        "β₁: $(round(F.beta1, sigdigits = 5)) ± $(round(F.beta1_se, sigdigits = 5))",
    )
    println(
        io,
        "x_intercept: $(round(F.x_intercept, sigdigits = 5)) ± $(round(F.x_intercept_se, sigdigits = 5))",
    )
    println(
        io,
        "χ²ᵣ (MSWD): $(round(F.reduced_chi_squared, sigdigits = 3)); p-value: $(round(F.p_value, sigdigits = 5))",
    )
    println(io, "no. observations: $(F.n_observations)")
    println(io, "FixedPoint: $(F.FixedPoint)")
end

struct GeneralisedLeastSquares <: LinearRegression
    beta::AbstractVector
    beta_se::AbstractVector
    variance_covariance::AbstractArray
    r_squared::AbstractFloat
    rmse::AbstractFloat
    n_observations::Integer
end

function Base.show(io::IOContext, F::GeneralisedLeastSquares)
    for i in eachindex(F.btea)
        println(io, "β$(i - 1): $(round.(F.beta[i]; digits = 5)) ± $(round.(F.beta_se[i]; digits = 5))")
    end
    println(io, "R²: $(round(F.r_squared; digits = 4))")
    println(io, "RMSE = $(round(F.rmse; digits=4))")
    return println(io, "no. observations: $(F.n_observations)")
end


function affine_fit(x::AbstractVector, fit::ErrorsInVariablesRegression)
    return fit.beta0 .+ x .* fit.beta1
end

function affine_ci(
    x::AbstractVector,
    fit::ErrorsInVariablesRegression;
    ci_level::AbstractFloat = 0.95,
)
    ν = fit.n_observations - 2
    tvalue = cquantile(TDist(ν), (1 - ci_level) / 2)
    return vec(sqrt.(abs.(fit.reduced_chi_squared  .* (x * fit.covariance_beta)))) .* tvalue
end

function affine_pi(
    x::AbstractVector,
    fit::ErrorsInVariablesRegression;
    ci_level::AbstractFloat = 0.95,
)
    ν = fit.n_observations - 2
    tvalue = cquantile(TDist(ν), (1 - ci_level) / 2)
    return vec(sqrt.(abs.(fit.reduced_chi_squared .* (1 .+ x .* fit.covariance_beta)))) .* tvalue
end
