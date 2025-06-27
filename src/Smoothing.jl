# SPDX-FileCopyrightText: 2025 Jarred Lloyd (https://github.com/jarredclloyd)
# SPDX-License-Identifier: MIT

# GeochemistryTools.jl/Smoothing.jl
# smoothing functions

export whittaker_smooth

"""
    whittaker_smooth(y::AbstractVector; kwargs)

Smooth some data `y` with equal spacing in `x`.

Input df as a DataFrame and specify the x and y data column names as symbols.

# Keywords

  - `lambda::Integer`: bandwidth parameter for smoothing.

        + Default value is UInt8(cld(sqrt(length(y)), 2))

  - `nthdiff::Integer`: Difference order to compute

        + Default value is `2` (second-order difference)

# Description

This function performs smoothing of `y` values that have equispaced `x` and are non-monotonic.

Algorithm uses the "Whittaker" smoother presented in Eilers (2003)

# References

Eilers, PHC (2003) 'A Perfect Smoother', *Analytical Chemistry*, 75(14):3631–3636. https://doi.org/10.1021/ac034173t
"""
function whittaker_smooth(
    y::AbstractVector;
    lambda::Integer = Integer(cld(sqrt(length(y)), 2)),
    nthdiff::Integer = 2,
)
    m = length(y)
    𝐈 = I(m)
    𝐃 = sparse(diff(𝐈; dims = 1))
    for i ∈ 2:nthdiff
        𝐃 = sparse(diff(𝐃; dims = 1))
    end
    Cm = cholesky(𝐈 + lambda * 𝐃' * 𝐃)
    return z = Cm \ (Cm' \ y)
end
