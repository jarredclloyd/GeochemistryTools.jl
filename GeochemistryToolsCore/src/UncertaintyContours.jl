# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

export uncertainty_ellipse, uncertainty_ellipsoid, mahalanobis_contour

"""
    mahalanobis_contour(μ, Σ; confidence_level=0.95, n=100, dims=nothing)

Return points on the Mahalanobis contour.

  - μ ∈ ℝⁿ
  - Σ ∈ ℝⁿˣⁿ
  - dims: optional index vector for projection (e.g. (1,2))
"""
function mahalanobis_contour(
    μ::AbstractVector{<:Real},
    Σ::AbstractMatrix{<:Real};
    confidence_level = 0.95,
    n = 100,
    dims
)
    dims = collect(dims)
    d = length(μ)
    @assert length(dims) == 2

    μ, A = mahalanobis_transform(μ, Σ; confidence_level)

    t = range(0, 2π; length = n)
    U = zeros(d, n)
    U[dims[1], :] .= cos.(t)
    U[dims[2], :] .= sin.(t)

    X = μ .+ A * U
    return X[dims, :]
end


function mahalanobis_transform(
    μ::AbstractVector{<:Real},
    Σ::AbstractMatrix{<:Real};
    confidence_level::Float64 = 0.95
)
    d = length(μ)
    @assert size(Σ) == (d, d)

    q = sqrt(quantile(Chisq(d), confidence_level))

    F = eigen(Symmetric(Σ))
    Σhalf = F.vectors * Diagonal(sqrt.(F.values)) * transpose(F.vectors)

    return μ, q * Σhalf
end

function covmat_rho(σ𝑥, σ𝑦, ρ𝑥𝑦; uncertainty_level_in = 2)
    σx = σ𝑥 / uncertainty_level_in
    σy = σ𝑦 / uncertainty_level_in
    σxy = ρ𝑥𝑦 * σx * σy
    return @SMatrix [
        σx^2  σxy
        σxy   σy^2
    ]
end

function uncertainty_ellipse(
    x::Measurement,
    y::Measurement;
    uncertainty_level_in::Integer = 2,
    confidence_level::Float64 = 0.95,
    n::Integer = 100
)
    return uncertainty_ellipse(
        Measurements.value(x),
        Measurements.value(y),
        Measurements.uncertainty(x),
        Measurements.uncertainty(y),
        Measurements.cor(x, y);
        uncertainty_level_in = uncertainty_level_in,
        confidence_level = confidence_level,
        n = n
    )
end

function uncertainty_ellipse(
    x::Real,
    y::Real,
    σx::Real,
    σy::Real,
    ρxy::Real;
    uncertainty_level_in::Integer = 2,
    confidence_level::Float64 = 0.95,
    n::Integer = 100
)
    μ = [x, y]
    Σ = covmat_rho(σx, σy, ρxy; uncertainty_level_in = uncertainty_level_in)

    return Matrix(transpose(mahalanobis_contour(
        μ,
        Σ;
        confidence_level = confidence_level,
        n = n,
        dims = (1, 2)
    )))
end


function uncertainty_ellipsoid(
    μ::AbstractVector{<:Measurement};
    confidence_level = 0.95,
    nθ = 50,
    nφ = 25
)
    return uncertainty_ellipsoid(
        Measurements.value.(μ),
        Measurements.cov(μ);
        confidence_level = confidence_level,
        nθ = nθ,
        nφ = nφ
    )
end

function uncertainty_ellipsoid(
    μ::AbstractVector{<:Real},
    Σ::AbstractMatrix{<:Real};
    confidence_level = 0.95,
    nθ = 50,
    nφ = 25
)
    d = length(μ)
    @assert d == 3

    μ, A = mahalanobis_transform(μ, Σ; confidence_level)

    θ = range(0, 2π; length = nθ)
    φ = range(0, π; length = nφ)

    X = zeros(3, nθ, nφ)
    R = zeros(nθ, nφ)

    for i ∈ eachindex(θ), j ∈ eachindex(φ)
        u = @SVector [cos(θ[i]) * sin(φ[j]), sin(θ[i]) * sin(φ[j]), cos(φ[j])]
        X[:, i, j] = μ .+ A * u

        Au = A * u
        R[i, j] = norm(Au)

    end

    return X, R   # 3 × nθ × nφ
end




# attempt at Makie recipe
# @recipe(ErrorEllipse) do scene
#     return Attributes(; confidence_level::Float64 = 0.95, n = 100)
# end


# function Makie.plot!(
#     p::ErrorEllipse,
#     centre::AbstractVector{<:Real},
#     Σ::AbstractMatrix{<:Real}
# )
#     pts = ellipse_points(centre, Σ; confidence_level = p.confidence_level[], n = p.n[])

#     return lines!(p, pts)
# end
