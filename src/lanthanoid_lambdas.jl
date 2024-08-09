#= Preamble
This file contains functions to load lanthanoid elemental data (+Y) and calculate lambda fits.
=#

function LambdaREE(
    lanth_measured::AbstractVector{AbstractString},
    lanth_values::AbstractVector{T} where T <: Real,
    lanth_uncertainties::Union{Nothing, AbstractArray} = nothing;
    weight_type::AbstractString = "abs",
    normalise::Bool = true,
    normalisation_values::AbstractString = "ci-chondrite, PO2016",
    fit_ce::Bool = false,
    fit_eu::Bool = false,
    fit_gd::Bool = true,
)
    lanthanoids =
        ["La", "Ce", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
    lanth_radii = deepcopy(cn_eight_IR.(lanthanoids .* "_3+"))
    lanth_uncertainties = lanth_uncertainties[isfinite.(lanth_values) .== true]
    lanth_values = lanth_values[isfinite.(lanth_values) .== true]
    𝑁::Integer = length(lanth_radii)
    if 𝑁 == length(lanth_values) && 𝑁 > 2
        x_sums::Vector{MultiFloat{Float64,4}} = Vector{MultiFloat{Float64,4}}(undef, 7)
        @simd for i ∈ eachindex(x_sums)
            x_sums[i] = sum(lanth_radii .^ i)
        end
        β::MultiFloat{Float64,4}         = _beta_orthogonal(𝑁, x_sums)
        γ::Vector{MultiFloat{Float64,4}} = _gamma_orthogonal(𝑁, x_sums)
        δ::Vector{MultiFloat{Float64,4}} = _delta_orthogonal(𝑁, x_sums)
        ϵ::Vector{MultiFloat{Float64,4}} = _epsilon_orthogonal(𝑁, x_sums)
        X::Matrix{MultiFloat{Float64,4}} = hcat(fill(1.0, 𝑁), (lanth_radii .- β), (lanth_radii .- γ[1]) .* (lanth_radii .- γ[2]), (lanth_radii .- δ[1]) .* (lanth_radii .- δ[2]) .* (lanth_radii .- δ[3]), (lanth_radii .- ϵ[1]) .* (lanth_radii .- ϵ[2]) .* (lanth_radii .- ϵ[3]) .* (lanth_radii .- ϵ[4]))
        X                                = X[1:size(X, 1) .!= findfirst(occursin("Pm"), lanthanoids), :]
        deleteat!(lanthanoids, occursin.("Pm", lanthanoids))

        measured_elements::Vector{AbstractString} = names(input_data)
        if normalise == true
            norm_values::Vector{Real} = Vector{Real}(undef, length(measured_elements))
            if lowercase(normalisation_values) == "ci-chondrite, PO2016"
                for i in eachindex(measured_elements)
                    norm_values[i] = deepcopy(ci_chondrite_PO2016(measured_elements[i])[1])
                end
                # elseif lowercase(normalisation_values) == "PAAS"
                #     for i in eachindex(measured_elements)
                #         norm_values[i] = deepcopy(PAAS(measured_elements[i])[1])
                #     end
            end
            normalised_data::Matrix{MultiFloat{Float64,4}} =
                transpose(input_data) ./ norm_values
        end

        if fit_gd == false
            X = X[1:size(X, 1) .!= findfirst(occursin("Gd"), lanthanoids), :]
        end
        if fit_eu == false
            X = X[1:size(X, 1) .!= findfirst(occursin("Eu"), lanthanoids), :]
        end
        if fit_ce == false
            X = X[1:size(X, 1) .!= findfirst(occursin("Ce"), lanthanoids), :]
        end
        order::Vector{Integer} = [0, 1, 2, 3, 4]
        if lanth_uncertainties === nothing
            ω::Vector{MultiFloat{Float64,4}} = fill(1.0, length(y))
        elseif occursin("rel", lowercase(weight_type)) === true
            ω = lanth_uncertainties
        elseif occursin("abs", lowercase(weight_type)) == true
            ω = abs.(lanth_uncertainties) ./ abs.(y)
        else
            throw(
                ArgumentError(
                    "Value of 'weight_type' is unrecognised. String should contain either 'rel' or 'abs'.",
                ),
            )
        end
        Ω::Diagonal{MultiFloat{Float64,4},Vector{MultiFloat{Float64,4}}} =
            Diagonal(1 ./ (ω ./ mean(ω)) .^ 2)
        Xᵀ::Transpose{MultiFloat{Float64,4},Matrix{MultiFloat{Float64,4}}} = transpose(X)
        VarΛX::Symmetric{Float64,Matrix{Float64}} = Symmetric(inv(Xᵀ * (Ω) * X))
        Λ::Vector{Float64} = VarΛX * Xᵀ * Ω * y
        for i in eachindex(Λ)
            Λ[i] = abs(Λ[i]) ≤ Base.rtoldefault(Float64) ? 0.0 : Λ[i]
        end
        rss::Vector{Float64} = Vector{Float64}(undef, 5)
        @inbounds @simd for i ∈ eachindex(order)
            residuals::Vector{MultiFloat{Float64,4}} = (y .- (view(X, :, 1:i) * Λ[1:i]))
            rss[i] = transpose(residuals) * Ω * (residuals)
        end
        mse = rss ./ (𝑁 .- (order .+ 1))
        Λ_SE::AbstractMatrix{Float64} = zeros(Float64, 5, 5)
        @inbounds for i ∈ eachindex(order)
            Λ_SE[1:i, i] = sqrt.(diag(view(VarΛX, 1:i, 1:i) * (mse[i])))
        end
        sparse(Λ_SE)
        tss::Float64 = transpose((y .- mean(y))) * Ω * (y .- mean(y))
        rmse::Vector{Float64} = sqrt.(mse)
        R²::Vector{Float64} = 1 .- (rss ./ (tss))
        R²ₒₚ::Vector{Float64} = _olkin_pratt.(R², 𝑁, order .+ 1)
        BIC::Vector{Float64} = Vector{Float64}(undef, 5)
        AIC::Vector{Float64} = Vector{Float64}(undef, 5)
        BIC = _bayesian_information_criteria.(rss, 𝑁, order)
        BICw =
            exp.(-0.5 .* (BIC .- minimum(BIC))) ./ sum(exp.(-0.5 .* (BIC .- minimum(BIC))))
        AIC = _akaike_information_criteria.(rss, 𝑁, order)
        AICw =
            exp.(-0.5 .* (AIC .- minimum(AIC))) ./ sum(exp.(-0.5 .* (AIC .- minimum(AIC))))
        return OrthogonalPolynomial(
            Λ,
            Λ_SE,
            Float64.(β),
            Float64.(γ),
            Float64.(δ),
            Float64.(ϵ),
            VarΛX,
            order,
            R²,
            R²ₒₚ,
            rmse,
            rss,
            mse,
            AIC,
            AICw,
            BIC,
            BICw,
            𝑁,
        )
    else
        println("Unable to fit data as there are less than three values")
        return OrthogonalPolynomial(
            fill(nothing, length(fieldnames(OrthogonalPolynomial)))...,
        )
    end
end
