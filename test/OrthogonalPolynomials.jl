using Test, Random, GeochemistryTools, Distributions, LinearAlgebra, MultiFloats, StatsBase, GenericLinearAlgebra
@testset "orthogonal polynomials"  begin
    # setup
    test_x = Float64x4.(collect(-1000000:1000:1000000));
    𝑁 = length(test_x)
    x_sums = Vector{MultiFloat{Float64,4}}(undef, 7)
            @simd for i ∈ eachindex(x_sums)
                x_sums[i] = sum(test_x.^ i)
            end
    β = GeochemistryTools._beta_orthogonal(𝑁, x_sums)
    γ = GeochemistryTools._gamma_orthogonal(𝑁, x_sums)
    δ = GeochemistryTools._delta_orthogonal(𝑁, x_sums)
    ϵ = GeochemistryTools._epsilon_orthogonal(𝑁, x_sums)
    order= [0, 1, 2, 3, 4]
    X = hcat(fill(1.0, 𝑁), (test_x .- β), (test_x .- γ[1]) .* (test_x .- γ[2]), (test_x .- δ[1]) .* (test_x .- δ[2]) .* (test_x .- δ[3]), (test_x .- ϵ[1]) .* (test_x .- ϵ[2]) .* (test_x .- ϵ[3]) .* (test_x .- ϵ[4]))

    per_err = abs.(rand(Xoshiro(), Normal(0.02,0.01), 𝑁));

    # testing
    Λ = [0.1, -0.1, 0.1, -0.1, 0.1];
    test_y = GeochemistryTools._poly_orthogonal(test_x, Λ, β, γ, δ, ϵ, 4);
    test_array = hcat(test_x, test_y);
    test_errors = per_err .* test_y;
    test_fit = fit_orthogonal(test_array);
    test_fit_wt = fit_orthogonal(hcat(test_array, test_errors); errors=true);
    @test isapprox(test_fit.lambda, Λ; rtol = sqrt(eps(Float64)))
    @test isapprox(test_fit_wt.lambda, Λ; rtol = sqrt(eps(Float64)))
    @test isapprox(poly_orthogonal(test_x, test_fit, 4), test_y; rtol = sqrt(eps(Float64)))

    Λ = collect(rand(Xoshiro(), Normal(0,20), 5));
    test_y = GeochemistryTools._poly_orthogonal(test_x, Λ, β, γ, δ, ϵ, 4);
    test_array = hcat(test_x, test_y);
    test_errors = per_err .* test_y;
    test_fit = fit_orthogonal(test_array);
    test_fit_wt = fit_orthogonal(hcat(test_array, test_errors); errors=true);
    @test isapprox(test_fit_wt.lambda, Λ; rtol = sqrt(eps(Float64)))
    @test isapprox(test_fit.lambda, Λ; rtol = sqrt(eps(Float64)))
    @test isapprox(poly_orthogonal(test_x, test_fit, 4), test_y; rtol = sqrt(eps(Float64)))
end
