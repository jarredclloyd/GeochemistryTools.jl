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
    println(io, "βₒ: $(F.beta0)")
    println(io, "βₒSE: $(F.beta0_se)")
    println(io, "beta1: $(F.beta1)")
    println(io, "beta1_se: $(F.beta1_se)")
    println(io, "x_intercept: $(F.x_intercept)")
    println(io, "x_intercept_se: $(F.x_intercept_se)")
    println(io, "reduced_chi_squared (MSWD): $(F.reduced_chi_squared)")
    return println(io, "no. observations: $(F.n_observations)")
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
    println(io, "βₒ: $(F.beta0)")
    println(io, "βₒSE: $(F.beta0_se)")
    println(io, "beta1: $(F.beta1)")
    println(io, "beta1_se: $(F.beta1_se)")
    println(io, "x_intercept: $(F.x_intercept)")
    println(io, "x_intercept_se: $(F.x_intercept_se)")
    println(io, "reduced_chi_squared (MSWD): $(F.reduced_chi_squared)")
    println(io, "no. observations: $(F.n_observations)")
    return println(io, "FixedPoint: $(F.FixedPoint)")
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
    println(io, "β₀: $(round.(F.beta; digits = 5)) ± $(round.(F.beta_se; digits = 5))")
    println(io, "R²: $(round(F.r_squared; digits = 4))")
    println(io, "RMSE = $(round(F.rmse; digits=4))")
    return println(io, "no. observations: $(F.n_observations)")
end


# function linear_fit(
#     df::AbstractDataFrame,
#     method::AbstractString = "Mahon";
#     [se_level_in::Integer = 2, se_level_out::Integer = 2, se_type::AbstractString = "abs", initial::Any = nothing],
# )

# end
