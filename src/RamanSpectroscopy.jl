#= Preamble
This file contains various decay system equations for calculation of ratios and ages (U-Pb).
=#

export anmirpls

#Functions
"""
    ANMIRPLS(ṽ, Ι)

Compute the baseline subtracted raman-shift spectra.

ṽ is a vector of the wavenumbers and Ι is a vector of intensities.
Uses the adaptive noise model based iteratively reweighted penalized least squares algorithm (ANM-IRPLS) of
Saveliev et al. 2022 (https://doi.org/10.1002/jrs.6275). 
"""
function anmirpls(ṽ, Ι)
    #1. basis construction
    fᵢ(x) = exp(-0.5 * (x - mᵢ)^2/σ^2)
    #2. regularisation

    #3. weights
    
    ωᵢ = 1/(10e-4 .+ Ys₁ .+ 10 .* Ys₂)^2
    #4. matrix for coefficient calculation
    
    #5. baseline approximation - iterative


    #6. correction of vertical displacement with noise
    m = median(Y₀ .- fit)
    bias = median((Y₀ .- fit) .* ind)
    #7. final baseline subtraction
    corrected = Y₀ .- fit - bias
end

"""
    IarPLS(v̄, I)

    Compute the baseline subtracted raman-shift spectra via IarPLS.

ṽ is a vector of the wavenumbers and Ι is a vector of intensities.
"""
function IarPLS(I, λ = 1e5, diff = 2, max_iter = 50, ϵ = 1e-6)
    n = length(I)
    t = 1
    Wⁱ = LinearAlgebra.I(n)
    for t ∈ 1:max_iter
        Z = (W + λ * D' * D) ^-1 * W * Y

end

function ISRU(x)
    1 / 2 * (1 - (ℯ^t * (x - 2 * σₓ₋) / σₓ₋))
    
end