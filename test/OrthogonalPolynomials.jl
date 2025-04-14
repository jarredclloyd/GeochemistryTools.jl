using Random, GeochemistryTools, Distributions, LinearAlgebra, MultiFloats, StatsBase
test_x = collect(0:0.5:60)
𝑁::Integer = length(test_x)

x_sums::Vector{MultiFloat{Float64,4}} = Vector{MultiFloat{Float64,4}}(undef, 7)
        @simd for i ∈ eachindex(x_sums)
            x_sums[i] = sum(test_x .^ i)
        end
β::MultiFloat{Float64,4} = GeochemistryTools._beta_orthogonal(𝑁, x_sums)
γ::Vector{MultiFloat{Float64,4}} = GeochemistryTools._gamma_orthogonal(𝑁, x_sums)
δ::Vector{MultiFloat{Float64,4}} = GeochemistryTools._delta_orthogonal(𝑁, x_sums)
ϵ::Vector{MultiFloat{Float64,4}} = GeochemistryTools._epsilon_orthogonal(𝑁, x_sums)
order::Vector{Integer} = [0, 1, 2, 3, 4]

X::Matrix{MultiFloat{Float64,4}} = hcat(fill(1.0, 𝑁), (test_x .- β), (test_x .- γ[1]) .* (test_x .- γ[2]), (test_x .- δ[1]) .* (test_x .- δ[2]) .* (test_x .- δ[3]), (test_x .- ϵ[1]) .* (test_x .- ϵ[2]) .* (test_x .- ϵ[3]) .* (test_x .- ϵ[4]))

Λ = [1, -1, 1, -1, 1]
test_y = GeochemistryTools._poly_orthogonal(test_x, Λ, β, γ, δ, ϵ, 4)
test_y_noisy = GeochemistryTools._poly_orthogonal(test_x .+ noise_x, Λ, β, γ, δ, ϵ, 4)
test_array = hcat(test_x, test_y)
test_errors = abs.(rand(Xoshiro(), Normal(0.02,0.01), 𝑁))
test_fit = fit_orthogonal(test_array)

@test test_fit.lambda .≈ Λ

noise_x = rand(Xoshiro(), Normal(0, 0.1),𝑁)
noise_y = rand(Xoshiro(), Normal(1, 0.1), 𝑁);
test_x_noisy = test_x .+ noise_x
test_y_noisy = test_y .* noise_y;
test_array_noisy = hcat(test_x_noisy, test_y_noisy);
test_fit_noisy = fit_orthogonal(test_array_noisy);
coeffs = reshape(test_fit_noisy.lambda, 1, 5);

for i ∈ 2:10000
    noise_x = rand(Xoshiro(), Normal(0, 0.1),𝑁)
    noise_y = rand(Xoshiro(), Normal(1, 0.1), 𝑁);
    test_x_noisy = test_x .+ noise_x
    test_y_noisy = test_y .* noise_y;
    test_array_noisy = hcat(test_x_noisy, test_y_noisy);
    test_fit_noisy = fit_orthogonal(test_array_noisy);
    coeffs = vcat(coeffs, test_fit_noisy.lambda')
    GLScoeffs = vcat(GLScoeffs, GLS_fit_noisy.beta')
end
coeff_mean_var = hcat(reshape(mean(coeffs; dims=1),5,1), reshape(var(coeffs; dims=1),5,1))
