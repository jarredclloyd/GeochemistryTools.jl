#= Preamble
This file contains various decay system equations for calculation of ratios and ages (U-Pb).
=#

export ageRbSr, RbSrAgeNorm, RbSrAgeInv

#Rb-Sr
"""
    ageRbSr(β₁, [β₀=0.0, β₁SE=0.0, β₀SE=0.0, σᵦ₁ᵦ₀=0.0; inverse = false, SElevel_in = 2, SElevel_out = 2])

    calculates ages and uncertainties of RbSr data.

"""
function ageRbSr(β₁, β₀=0.0, β₁SE=0.0, β₀SE=0.0, σᵦ₁ᵦ₀=0.0; 
    inverse = false, SElevel_in = 2, SElevel_out = 2)
    if inverse == false
        if .==(β₁, 0.0) == true
            error("Non-zero values ∈ ℝ are required for slope (β₁) to calculate ages for isochron data.")
        end
        if isa(β₁, Array) == true
            nslopes = length(β₁)
            age = zeros(nslopes)
            ageSE = zeros(nslopes)
            if isa(β₁SE, Array) == false
                β₁SE = zeros(nslopes)
            end
            @simd for i in 1:nslopes
                age[i], ageSE[i] = RbSrAgeNorm(β₁[i], β₁SE[i] / SElevel_in)
                ageSE[i] = ageSE[i] * SElevel_out
            end
        else 
            age, ageSE = RbSrAgeNorm(β₁, β₁SE / SElevel_in)
            ageSE = ageSE * SElevel_out
        end
    elseif inverse == true
        if .==(β₀, 0.0) == true || .==(β₁, 0.0) == true
            error("Non-zero values ∈ ℝ are required for both the intercept (β₀) and the slope (β₁) to calculate ages for 
            inverse isochron data.")
        end
        if isa(β₁, Array) == true && isa(β₀, Array) == true && length(β₁) == length(β₀)
            nslopes = length(β₁)
            age = zeros(nslopes)
            ageSE = zeros(nslopes)
            if isa(β₀SE, Array) == false
                β₀SE = zeros(nslopes)
            end
            if isa(β₁SE, Array) == false
                β₁SE = zeros(nslopes)
            end
            if isa(σᵦ₁ᵦ₀, Array) == false
                σᵦ₁ᵦ₀ = zeros(nslopes)
            end
            @simd for i in 1:nslopes
                age[i], ageSE[i] = RbSrAgeInv(β₀[i], β₁[i], β₀SE[i] / SElevel_in, β₁SE[i] / SElevel_in, σᵦ₁ᵦ₀[i])
                ageSE[i] = ageSE[i] * SElevel_out
            end
        else 
            age, ageSE = RbSrAgeInv(β₀, β₁, β₀SE / SElevel_in, β₁SE / SElevel_in)
            ageSE = ageSE * SElevel_out
        end
    end
    return age, ageSE
end 

function RbSrAgeNorm(β₁::Real, β₁SE::Real=0.0)
    if .==(β₁, 0.0) == true
        error("Non-zero values ∈ ℝ are required for slope (β₁) to calculate ages for isochron data.")
    end
    date = log(β₁ + 1) / λRb87
    dateSE = abs(log(β₁SE + 1) / λRb87)
    return date, dateSE
end

function RbSrAgeInv(β₀::Real, β₁::Real, β₀SE::Real=0.0, β₁SE::Real=0.0, σᵦ₁ᵦ₀::Real=0.0)
    if .==(β₀, 0.0) == true || .==(β₁, 0.0) == true
            error("Non-zero values ∈ ℝ are required for both the intercept (β₀) and the slope (β₁) to calculate ages for 
            inverse isochron data.")
    end
    ratio = inv(-β₀ / β₁)
    date = log(ratio + 1) / λRb87
    ratioSE = ratio * sqrt((β₀SE / β₀)^2 + (β₁SE / β₁)^2 - 2 * σᵦ₁ᵦ₀ / (β₀ * β₁))
    dateSE = abs(log(ratioSE + (β₀ ^2 * β₀SE ^2 + β₁ ^2 * β₁SE ^2 + 2 * β₀ * β₁ * σᵦ₁ᵦ₀) + 1) / λRb87)
    return date, dateSE
end