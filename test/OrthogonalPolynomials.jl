using Random, GeochemistryTools, Distributions
Random.seed!(12)
test_x = collect(0:0.5:30)
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

Λ = [10, 1, 0.1, 0.01, 0.001]
Λ2 = collect(rand(Xoshiro(12), 5))
test_ϵ = rand(Normal(0,0.01), 𝑁)

test_y = GeochemistryTools._poly_orthogonal(test_x, Λ, β, γ, δ, ϵ, 4)
test_y_noisy = GeochemistryTools._poly_orthogonal(test_x, Λ, β, γ, δ, ϵ, 4) .* (1 .+ test_ϵ)
test_y2 = GeochemistryTools._poly_orthogonal(test_x, Λ2, β, γ, δ, ϵ, 4) .* (1 .+ test_ϵ)


test_array = hcat(test_x, test_y)
test_array_noisy = hcat(test_x, test_y_noisy)
test_array2 = hcat(test_x, test_y2)

test_fit = fit_orthogonal(test_array)
test_fit_noisy = fit_orthogonal(test_array_noisy)
test_fit2 = fit_orthogonal(test_array2)

test_fit.lambda .≈ Λ

ω = abs.(test_ϵ)
Ω::Diagonal{MultiFloat{Float64,4},Vector{MultiFloat{Float64,4}}} =
Diagonal((ω ./ mean(ω)) .^ 2)
VarΛX::Symmetric{Float64,Matrix{Float64}} = Symmetric(inv(transpose(X) * inv(Ω) * X))
VarΛX * transpose(X) * inv(Ω) * y_noise
