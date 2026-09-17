# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

export olkin_pratt,
    chi_squared_reduced,
    akaike_information_criteria,
    bayesian_information_criteria,
    chi_squared_reduced_confint


"""
    olkin_pratt(R²::Real, 𝑛::Integer, predictors::Integer)

Compute the unbiased Olkin-Pratt estimator of ρ² from R².

```
julia> olkin_pratt(0.95, 100, 3)
0.9494275113877586
```

# References

Olkin, I, & Pratt, JW (1958) 'Unbiased Estimation of Certain Correlation Coefficients', *The Annals of Mathematical Statistics*, 29(1):201–211, https://doi.org/10.1214/aoms/1177706717

Karch, J (2020) 'Improving on Adjusted R-Squared', *Collabra: Psychology*, 6(1):45, https://doi.org/10.1525/collabra.343
"""
function olkin_pratt(R²::Real, 𝑛::Integer, predictors::Integer)
    if isfinite(R²) == false || R² < 0 || R² > 1
        return NaN
    else
        z = 1 - R²
        c = (𝑛 - predictors + 1) / 2
        if c ≤ 2
            return NaN
        else
            if z == 0
                _₂F₁value = 1
            elseif z == 1
                _₂F₁value = (c - 1) / (c - 2)
            else
                _₂F₁value = _2F1_taylor(1, 1, c, z)
            end
            return 1 - ((𝑛 - 3) / (𝑛 - predictors - 1)) * z * _₂F₁value
        end
    end
end

function chi_squared_reduced(χ²::Real, 𝑛::Integer, predictors::Integer)
    return χ² / (𝑛 - predictors)
end

function bayesian_information_criteria(rss::Real, 𝑛::Integer, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be ≥ 0"))
    end
    if rss ≤ 0 || !isfinite(rss)
        throw(ArgumentError("RSS must be positive and finite"))
    end
    if 𝑛 ≤ 0
        throw(ArgumentError("Sample size 𝑛 must be positive for BIC"))
    end
    𝑘 = order + 2
    return 𝑛 * log(rss / 𝑛) + 𝑘 * log(𝑛)
end


function akaike_information_criteria(rss::Real, 𝑛::Integer, order::Integer)
    if order < 0
        throw(ArgumentError("Polynomial order must be ≥ 0"))
    end
    if rss ≤ 0 || !isfinite(rss)
        throw(ArgumentError("RSS must be positive and finite"))
    end
    𝑘 = order + 2
    if 𝑛 ≤ 𝑘 + 1
        throw(ArgumentError("Sample size 𝑛 ($𝑛) must be > 𝑘 + 1 ($(𝑘 + 1)) for AIC"))
    end
    return 𝑛 * log(rss / 𝑛) + 2 * 𝑘 + ((2 * 𝑘 * (𝑘 + 1)) / (𝑛 - 𝑘 - 1))
end

function chi_squared_reduced_confint(dof::Integer, confidence_level::AbstractFloat = 0.95)
    lower_χ²ᵣ = cquantile(Chisq(dof), 1 - (1 - confidence_level) / 2) / dof
    upper_χ²ᵣ = cquantile(Chisq(dof), (1 - confidence_level) / 2) / dof
    return (lower_χ²ᵣ, upper_χ²ᵣ)
end

function _2F1_taylor(
    a::Real,
    b::Real,
    c::Real,
    z::Real;
    tol = eps(eltype(z)),
    maxiter::Int64 = 1000
)
    @assert a == b == 1 "This implementation of the ₂F₁ function is only valid for a == b == 1"
    T = eltype(z)
    Cⱼ, Sⱼ = T(1), T(1)
    j = 0
    while abs(Cⱼ) ≥ tol * (one(T) + abs(Sⱼ)) && j ≤ maxiter
        Cⱼ *= (j + 1) / (c + j) * z # = (a + j) * (b + j) / (c + j) * z / (j + 1), a=b=1
        Sⱼ += Cⱼ
        j += 1
    end
    return Float64(Sⱼ)
end
