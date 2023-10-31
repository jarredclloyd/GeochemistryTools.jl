#= Preamble
This file contains various decay system equations for calculation of ratios and ages (Rb—Sr).
=#

export ageRbSr, RbSrAgeNorm, RbSrAgeInv, confidence_interval, inflate_ci, inflate_se

#Rb-Sr
"""
    ageRbSr(β₁, [β₀=0.0, β₁_se=0.0, β₀_se=0.0, σᵦ₁ᵦ₀=0.0; inverse = false, se_level_in = 2, se_level_out = 2])

    calculates ages and uncertainties of RbSr data.

    # Arguments
"""
function ageRbSr(
    β₁::AbstractFloat,
    β₀::AbstractFloat = 0.0,
    β₁_se::AbstractFloat = 0.0,
    β₀_se::AbstractFloat = 0.0,
    σᵦ₁ᵦ₀::AbstractFloat = 0.0;
    inverse::Bool = false,
    se_level_in::Int = 1,
    se_level_out::Int = 2,
)
    if inverse == false
        if .==(β₁, 0.0) == true
            throw(
                ArgumentError(
                    "Non-zero values ∈ ℝ are required for slope (β₁) to calculate ages for isochron data.",
                ),
            )
        end
        if isa(β₁, Array) == true
            nslopes = length(β₁)
            age = zeros(nslopes)
            age_se = zeros(nslopes)
            if isa(β₁_se, Array) == false
                β₁_se = zeros(nslopes)
            end
            @simd for i = 1:nslopes
                age[i], age_se[i] = RbSrAgeNorm(β₁[i], β₁_se[i] / se_level_in)
                age_se[i] = age_se[i] * se_level_out
            end
        else
            age, age_se = RbSrAgeNorm(β₁, β₁_se / se_level_in)
            age_se = age_se * se_level_out
        end
    elseif inverse == true
        if .==(β₀, 0.0) == true || .==(β₁, 0.0) == true
            throw(
                ArgumentError(
                    "Non-zero values ∈ ℝ are required for both the intercept (β₀) and the slope (β₁) to calculate ages for inverse isochron data.",
                ),
            )
        end
        if isa(β₁, Array) == true && isa(β₀, Array) == true && length(β₁) == length(β₀)
            nslopes = length(β₁)
            age = zeros(nslopes)
            age_se = zeros(nslopes)
            if isa(β₀_se, Array) == false
                β₀_se = zeros(nslopes)
            end
            if isa(β₁_se, Array) == false
                β₁_se = zeros(nslopes)
            end
            if isa(σᵦ₁ᵦ₀, Array) == false
                σᵦ₁ᵦ₀ = zeros(nslopes)
            end
            @simd for i = 1:nslopes
                age[i], age_se[i] = RbSrAgeInv(
                    β₀[i],
                    β₁[i],
                    β₀_se[i] / se_level_in,
                    β₁_se[i] / se_level_in,
                    σᵦ₁ᵦ₀[i],
                )
                age_se[i] = age_se[i] * se_level_out
            end
        else
            age, age_se = RbSrAgeInv(β₀, β₁, β₀_se / se_level_in, β₁_se / se_level_in)
            age_se = age_se * se_level_out
        end
    end
    return age, age_se
end

function ageRbSr(fit::ErrorsInVariablesRegression; inverse = true, se_level_out::Integer = 2)
    if inverse == false
        age, age_se = RbSrAgeNorm(fit.beta1, fit.beta1_se)
        age_se = se_level_out * age_se
    else
        age, age_se = RbSrAgeNorm(inv(fit.x_intercept), inv(fit.x_intercept / fit.x_intercept_se))
        age_se = se_level_out * age_se
    end
    return age, age_se
end

function RbSrAgeNorm(β₁::AbstractFloat, β₁_se::AbstractFloat = 0.0)
    if .==(β₁, 0.0) == true
        throw(
            ArgumentError(
                "Non-zero values ∈ ℝ are required for slope (β₁) to calculate ages for isochron data.",
            ),
        )
    end
    date = log(β₁ + 1) / λRb87
    date_se = abs(log(β₁_se + 1) / λRb87)
    return date, date_se
end

function RbSrAgeInv(
    β₀::AbstractFloat,
    β₁::AbstractFloat,
    β₀_se::AbstractFloat = 0.0,
    β₁_se::AbstractFloat = 0.0,
    σᵦ₁ᵦ₀::AbstractFloat = 0.0,
)
    if .==(β₀, 0.0) == true || .==(β₁, 0.0) == true
        throw(
            ArgumentError(
                "Non-zero values ∈ ℝ are required for both the intercept (β₀) and the slope (β₁) to calculate ages for inverse isochron data.",
            ),
        )
    end
    ratio = inv(-β₀ / β₁)
    date = log(ratio + 1) / λRb87
    ratio_se = ratio * sqrt((β₀_se / β₀)^2 + (β₁_se / β₁)^2 - 2 * σᵦ₁ᵦ₀ / (β₀ * β₁))
    date_se = abs(
        log(ratio_se + (β₀^2 * β₀_se^2 + β₁^2 * β₁_se^2 + 2 * β₀ * β₁ * σᵦ₁ᵦ₀) + 1) / λRb87,
    )
    return date, date_se
end

function confidence_interval(se, n; se_level = 1, ci_level = 0.95)
    return (se / se_level) * cquantile(TDist(n - 2), (1 - ci_level) / 2)
end

function inflate_se(se, χ²ᵣ)
    return se * sqrt(χ²ᵣ)
end
function inflate_ci(se, χ²ᵣ, n; se_level = 1, ci_level = 0.95)
    return (se / se_level) * cquantile(TDist(n - 2), (1 - ci_level) / 2) * sqrt(χ²ᵣ)
end
