abstract type LinearRegression <:Any end
abstract type ErrorsInVariablesRegression <: LinearRegression end

Base.getindex(LR::LinearRegression, i::Integer) = getfield(LR, i)
Base.keys(LR::LinearRegression) = Base.keys(LR)

struct MahonNonFixed <: ErrorsInVariablesRegression
    β₀::Real
    β₀SE::Real
    β₁::Real
    β₁SE::Real
    X_intercept::Real
    X_interceptSE::Real
    χ²ᵣ::Real
    pval::Real
    N::Integer
end

function Base.show(io::IOContext, F::MahonNonFixed)
    println(io, "βₒ: $(F.β₀)")
    println(io, "βₒSE: $(F.β₀SE)")
    println(io, "β₁: $(F.β₁)")
    println(io, "β₁SE: $(F.β₁SE)")
    println(io, "X_intercept: $(F.X_intercept)")
    println(io, "X_interceptSE: $(F.X_interceptSE)")
    println(io, "χ²ᵣ (MSWD): $(F.χ²ᵣ)")
    return println(io, "no. observations: $(F.N)")
end

struct MahonFixed <: ErrorsInVariablesRegression
    β₀::Real
    β₀SE::Real
    β₁::Real
    β₁SE::Real
    X_intercept::Real
    X_interceptSE::Real
    χ²ᵣ::Real
    pval::Real
    N::Integer
    FixedPoint::Tuple{Real,Real,Real,Real}
end

function Base.show(io::IOContext, F::MahonFixed)
    println(io, "βₒ: $(F.β₀)")
    println(io, "βₒSE: $(F.β₀SE)")
    println(io, "β₁: $(F.β₁)")
    println(io, "β₁SE: $(F.β₁SE)")
    println(io, "X_intercept: $(F.X_intercept)")
    println(io, "X_interceptSE: $(F.X_interceptSE)")
    println(io, "χ²ᵣ (MSWD): $(F.χ²ᵣ)")
    println(io, "no. observations: $(F.N)")
    return println(io, "FixedPoint: $(F.FixedPoint)")
end

struct LeastSquares <: LinearRegression
    β₀::Real
    β₁::Real
    N::Integer
end

function Base.show(io::IOContext, F::LeastSquares)
    println(io, "βₒ: $(F.β₀)")
    println(io, "β₁: $(F.β₁)")
    println(io, "no. observations: $(F.N)")
end

function _least_squares(X::AbstractVector{T}, Y::AbstractVector{T}) where T<:Real
    if length(X) != length(Y)
        throw(DimensionMismatch("X and Y vectors are"))
    end
    N = length(X)
    β₁ = (N * sum(X .* Y) - sum(X) * sum(Y)) / (N * sum(X.^2) - (sum(X)^2))
    β₀ = (sum(Y) - β₁ * sum(X)) / N
    return LeastSquares(β₀, β₁, N)
end

# function linear_fit(
#     df::AbstractDataFrame,
#     method::AbstractString = "Mahon";
#     [se_level_in::Integer = 2, se_level_out::Integer = 2, se_type::AbstractString = "abs", initial::Any = nothing],
# )

# end
