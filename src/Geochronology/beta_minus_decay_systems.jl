#= Preamble
This file contains various decay system equations for calculation of ratios and ages (Rb—Sr).
=#

export ageRbSr, RbSrAgeNorm, RbSrAgeInv, confidenceInterval, inflateCI, inflateSE

#Rb-Sr
"""
    ageRbSr(β₁, [β₀=0.0, β₁SE=0.0, β₀SE=0.0, σᵦ₁ᵦ₀=0.0; inverse = false, se_level_in = 2, se_level_out = 2])

    calculates ages and uncertainties of RbSr data.

    # Arguments
"""
function ageRbSr(
    β₁::AbstractFloat,
    β₀::AbstractFloat = 0.0,
    β₁SE::AbstractFloat = 0.0,
    β₀SE::AbstractFloat = 0.0,
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
            if isa(β₁SE, Array) == false
                β₁SE = zeros(nslopes)
            end
            @simd for i = 1:nslopes
                age[i], age_se[i] = RbSrAgeNorm(β₁[i], β₁SE[i] / se_level_in)
                age_se[i] = age_se[i] * se_level_out
            end
        else
            age, age_se = RbSrAgeNorm(β₁, β₁SE / se_level_in)
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
            if isa(β₀SE, Array) == false
                β₀SE = zeros(nslopes)
            end
            if isa(β₁SE, Array) == false
                β₁SE = zeros(nslopes)
            end
            if isa(σᵦ₁ᵦ₀, Array) == false
                σᵦ₁ᵦ₀ = zeros(nslopes)
            end
            @simd for i = 1:nslopes
                age[i], age_se[i] = RbSrAgeInv(
                    β₀[i],
                    β₁[i],
                    β₀SE[i] / se_level_in,
                    β₁SE[i] / se_level_in,
                    σᵦ₁ᵦ₀[i],
                )
                age_se[i] = age_se[i] * se_level_out
            end
        else
            age, age_se = RbSrAgeInv(β₀, β₁, β₀SE / se_level_in, β₁SE / se_level_in)
            age_se = age_se * se_level_out
        end
    end
    return age, age_se
end

function ageRbSr(fit::ErrorsInVariablesRegression; inverse = true, se_level_out::Int = 2)
    if inverse == false
        age, age_se = RbSrAgeNorm(fit.beta1, fit.beta1_se)
        age_se = se_level_out * age_se
    else
        age, age_se = RbSrAgeNorm(inv(fit.x_intercept), inv(fit.x_intercept / fit.x_intercept_se))
        age_se = se_level_out * age_se
    end
    return age, age_se
end

function RbSrAgeNorm(β₁::AbstractFloat, β₁SE::AbstractFloat = 0.0)
    if .==(β₁, 0.0) == true
        throw(
            ArgumentError(
                "Non-zero values ∈ ℝ are required for slope (β₁) to calculate ages for isochron data.",
            ),
        )
    end
    date = log(β₁ + 1) / λRb87
    dateSE = abs(log(β₁SE + 1) / λRb87)
    return date, dateSE
end

function RbSrAgeInv(
    β₀::AbstractFloat,
    β₁::AbstractFloat,
    β₀SE::AbstractFloat = 0.0,
    β₁SE::AbstractFloat = 0.0,
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
    ratioSE = ratio * sqrt((β₀SE / β₀)^2 + (β₁SE / β₁)^2 - 2 * σᵦ₁ᵦ₀ / (β₀ * β₁))
    dateSE = abs(
        log(ratioSE + (β₀^2 * β₀SE^2 + β₁^2 * β₁SE^2 + 2 * β₀ * β₁ * σᵦ₁ᵦ₀) + 1) / λRb87,
    )
    return date, dateSE
end

function confidenceInterval(SE, n; se_level = 1, CIlevel = 0.95)
    return (SE / se_level) * cquantile(TDist(n - 2), (1 - CIlevel) / 2)
end

function inflateSE(SE, χ²ᵣ)
    return SE * sqrt(χ²ᵣ)
end
function inflateCI(SE, χ²ᵣ, n; se_level = 1, CIlevel = 0.95)
    return (SE / se_level) * cquantile(TDist(n - 2), (1 - CIlevel) / 2) * sqrt(χ²ᵣ)
end
