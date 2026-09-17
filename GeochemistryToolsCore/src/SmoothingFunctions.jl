# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

# GeochemistryToolsCore.jl/SmoothingFunctions.jl
# Provides a Whittaker smoothing function and a despiking algorithm.

export whittaker, whittaker!, despike, despike!

"""
    whittaker(y::AbstractVector{Real}; kwargs)

Smooth some data `y` with equal spacing in `x`.

# Keywords

  - `lambda::Real`: bandwidth parameter for smoothing.

      + Higher values will result in more aggressive smoothing

          * Default value is 1e2

  - `nthdiff::Integer`: Difference order to compute

      + Default value is `2` (second-order difference)

# Description

This function performs smoothing of `y` values that have equispaced `x` and are non-monotonic.

Algorithm uses the "Whittaker" smoother presented in Eilers (2003)

# References

Eilers, PHC (2003) 'A Perfect Smoother', *Analytical Chemistry*, 75(14):3631–3636. https://doi.org/10.1021/ac034173t
"""
function whittaker(y::AbstractVector{<:Real}; lambda::Real = 1e2, nthdiff::Integer = 2)
    m = length(y)
    𝐈 = I(m)
    D1(m) = spdiagm(m - 1, m, 0 => fill(-1.0, m - 1), 1 => fill(1.0, m - 1))
    𝐃 = D1(m)
    for _ ∈ 2:nthdiff
        𝐃 = D1(size(𝐃, 1)) * 𝐃
    end
    𝐂 = cholesky(𝐈 + lambda * 𝐃' * 𝐃)
    return 𝐂.U \ (𝐂.U' \ y)
end


"""
    whittaker!(y::AbstractVector{Real}; kwargs)

Smooth some data `y` with equal spacing in `x` in place.

**Warning: use with care as this will overwrite original data.**

Use `whittaker` to allocate a new vector with the result instead.
"""
function whittaker!(y::AbstractVector{<:Real}; lambda::Real = 1e2, nthdiff::Integer = 2)
    m = length(y)
    𝐈 = I(m)
    D1(m) = spdiagm(m - 1, m, 0 => fill(-1.0, m - 1), 1 => fill(1.0, m - 1))
    𝐃 = D1(m)
    for _ ∈ 2:nthdiff
        𝐃 = D1(size(𝐃, 1)) * 𝐃
    end
    𝐂 = cholesky(𝐈 + lambda * 𝐃' * 𝐃)
    return y .= 𝐂 \ (𝐂' \ y)
end

"""
    despike(y::AbstractVector{<:Real}; threshold::Real=6.0, bandwidth::Integer=5, interpolate::Bool=true)

Despike data in vector `y`.

# Description

Remove spikes in time-series data. Implements a simple algorithm outlined in
Whitaker & Hayes (2018) that uses modified z-scores to detect outlier points.

The argument `interpolate` can be set to false to return a vector of spike indices instead.

# References

Whitaker, DA & Hayes, K (2018) `A simple algorithm for despiking Raman spectra`,  *Chemometrics and Intelligent Laboratory Systems*, 179:82–84, https://doi.org/10.1016/j.chemolab.2018.06.009
"""
function despike(
    y::AbstractVector{<:Real};
    threshold::Real = 6.0,
    bandwidth::Integer = 5,
    interpolate::Bool = true
)
    despiked_y = copy(y)
    z_score::Vector{Float64} = diff(y)
    if length(z_score) < 2
        return interpolate ? despiked_y : Int[]
    end

    zy_median::Real = median(z_score)

    zy_mad::Real = mad(z_score; normalize = true) # normalise accounts for bias in MAD calculation

    z_score .= abs.(0.6745 .* (z_score .- zy_median) ./ zy_mad)

    spikes::Vector{Int64} = [findall(>(threshold), z_score) .+ 1...]

    if interpolate === true
        for i ∈ spikes
            despiked_y[i] = mean(y[[
                collect(max(firstindex(z_score), i - bandwidth):(i - 1))...,
                collect((i + 1):min(lastindex(z_score), i + bandwidth))...
            ]],)
        end
        return despiked_y
    else
        return spikes
    end
end


"""
    despike!(y::AbstractVector{<:Real}; threshold::Real=6.0, bandwidth::Integer=5)

Despike data in vector `y` in place.

**Warning: use with care as this will overwrite original data.**

Use `despike` to allocate a new vector with the result instead.

# References

Whitaker, DA & Hayes, K (2018) `A simple algorithm for despiking Raman spectra`,  *Chemometrics and Intelligent Laboratory Systems*, 179:82–84, https://doi.org/10.1016/j.chemolab.2018.06.009
"""
function despike!(y::AbstractVector{<:Real}; threshold::Real = 6.0, bandwidth::Integer = 5)

    z_score::Vector{Float64} = diff(y)
    if length(z_score) < 2
        return interpolate ? y : Int[]
    end

    zy_median::Real = median(z_score)

    zy_mad::Real = mad(z_score; normalize = true) # normalise accounts for bias in MAD calculation

    z_score .= abs.(0.6745 .* (z_score .- zy_median) ./ zy_mad)

    spikes::Vector{Int64} = [findall(>(threshold), z_score) .+ 1...]

    if interpolate === true
        for i ∈ spikes
            y[i] = mean(y[[
                collect(max(firstindex(z_score), i - bandwidth):(i - 1))...,
                collect((i + 1):min(lastindex(z_score), i + bandwidth))...
            ]],)
        end
        return y
    else
        return spikes
    end
end
