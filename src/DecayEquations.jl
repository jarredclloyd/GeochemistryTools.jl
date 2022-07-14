#=
This file contains various decay system equations for calculation of ratios and ages (U-Pb).
=#

export ratioPb207U235, ratioPb206U238, ratioPb207Pb206, agePb207U235, agePb206U238, agePb207Pb206

# U-Pb decay system
function ratioPb207U235(age)
    exp(lambdaU235 * age) - 1
end

function ratioPb206U238(age)
    exp(lambdaU238 * age) - 1
end

function ratioPb207Pb206(age)
    ((exp(lambdaU235 * age) - 1) / (exp(lambdaU238 * age) - 1)) * inv(137.818)
end

function ratioprimePb207Pb206(age) #  d/dt(((e^(lambdaU235 * age) - 1)/(e^(lambdaU238 * age) - 1))/137.818)
    inv(ratioU238U235) * ((lambdaU235) * exp(lambdaU235 * age) * (exp(lambdaU238 * age) - 1) - lambdaU238 * 
    exp(lambdaU238 * age) * (exp(lambdaU235 * age) - 1)) / (exp(lambdaU238 * age) - 1)^2
end

function agePb206U238(ratio)
    log(ratio + 1) / lambdaU238
end

function agePb207U235(ratio)
    log(ratio + 1) / lambdaU235
end

function agePb207Pb206(ratio)
    t_guess = (log(ratio) + 3.21121) / 0.000586671
    tolerance = 0.000000000001
    ε = 1*10^-24
    while abs(ratio - ratioPb207Pb206(t_guess)) > tolerance && abs(ratioprimePb207Pb206(t_guess)) > ε
        t_iteration = t_guess - (ratioPb207Pb206(t_guess) - ratio) / ratioprimePb207Pb206(t_guess)
        t_guess = t_iteration
    end
    return t_guess
end