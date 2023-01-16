#= Preamble
This file contains various decay system equations for calculation of ratios and ages (U-Pb).
=#

export ageRbSr, RbSrAgeNorm, RbSrAgeInv

#Rb-Sr

function ageRbSr(β₁::Real, β₁SE::Real=0, β₀::Real=0, β₀SE::Real=0, σᵦ₁ᵦ₀::Real=0; inverse = false, SElevel_in = 2, 
    SElevel_out = 2)
    if inverse == false
        if isa(β₁, Array) == true
        age = zeros(size(β1))
        ageSE = zeros(size(β1SE))
        nslopes = length(β₁)
        @simd for i in 1:nslopes
            age[i], ageSE[i] * SElevel_out = RbSrAgeNorm(β₁[i], β₁SE[i] / SElevel_in)
            end
        else 
            age, ageSE * SElevel_out = RbSrAgeNorm(β₁[i], β₁SE[i] / SElevel_in)
        end
    elseif inverse == true
        if isa.(β₀, Real) !== true || if isa.(β₁, Real) !== true
            error("Both the intercept (β₀) and the slope (β₁) are required for inverse isochron age calculations.")
        end
        if isa(β₁, Array) == true && isa(β₀, Array) == true && length(β₁) == length(β₀)
            age = zeros(size(β1))
            ageSE = zeros(size(β1SE))
            nslopes = length(β₁)
            @simd for i in 1:nslopes
                age[i], ageSE[i] * SElevel_out = RbSrAgeInv(β₀[i], β₁[i], β₀SE[i] / SElevel_in, β₁SE[i] / SElevel_in, σᵦ₁ᵦ₀[i])
                end
        else 
            age, ageSE * SElevel_out = RbSrAgeInv(β₀[i], β₁[i], β₀SE[i] / SElevel_in, β₁SE[i] / SElevel_in)
        end
    end
    return age, ageSE
end

function RbSrAgeNorm(β₁::Real, β₁SE::Real=0)
    date = log(β₁ + 1) / λRb87
    dateSE = abs(log(β₁SE + 1) / λRb87)
    if dateSE > 0
        return date, dateSE
    else
        return date
    end
end

function RbSrAgeInv(β₀::Real, β₁::Real, β₀SE::Real=0, β₁SE::Real=0, σᵦ₁ᵦ₀::Real=0)
    if isa(β₀, Real) !== true || if isa(β₁, Real) !== true
        error("Both the intercept (β₀) and the slope (β₁) are required for inverse isochron age calculations.")
    end
    ratio = inv(-β₀ / β₁)
    date = log(ratio + 1) / λRb87
    ratioSE = ratio * sqrt((β₀SE / β₀)^2 + (β₁SE / β₁)^2 - 2 * σᵦ₁ᵦ₀ / (β₀ * β₁))
    dateSE = abs(log(ratioSE + (β₀ ^2 * β₀SE ^2 + β₁ ^2 * β₁SE ^2 + 2 * β₀ * β₁ * σᵦ₁ᵦ₀) + 1) / λRb87)
    if dateSE > 0
        return date, dateSE
    else
        return date
    end
end