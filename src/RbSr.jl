#= Preamble
This file contains various decay system equations for calculation of ratios and ages (U-Pb).
=#

export RbSrAgeNorm, RbSrAgeInv

#Rb-Sr

function RbSrAgeNorm(β₁)
    date = log(β₁+1)/λRb87    
end

function RbSrAgeInv(β₁)
    date = log(inv(β₁)+1)/λRb87   
end