#= Preamble
This file contains functions to load lanthanoid elemental data (+Y) and calculate lambda fits.
=#

function LambdaREE(
    lanth_measured::AbstractVector{AbstractString},
    lanth_values::AbstractVector{T} where T <: Real,
    lanth_uncertainties::Union{Nothing, AbstractArray} = nothing;
    weight_type::AbstractString = "abs",
    normalise::Bool = true,
    normalisation_values::AbstractString = "CI_PO2016",
    fit_ce::Bool = false,
    fit_eu::Bool = false,
    fit_gd::Bool = true,
)
    lanthanoids =
        ["La", "Ce", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
    lanth_radii = IONIC_RADIUS_CNEIGHT.(lanthanoids .* "_3+")
    if isnothing(lanth_uncertainties) != true
        finite_indices = findall(isfinite, lanth_values)
        # lanth_values =
    else
        finite_indices = intersect(findall(isfinite, lanth_values), findall(isfinite, lanth_measured))
    end


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

        X = X[1:size(X, 1) .!= findfirst(occursin("Pm"), lanthanoids), :]
        deleteat!(lanthanoids, occursin.("Pm", lanthanoids))

        measured_elements::Vector{AbstractString} = names(input_data)
        if normalise == true
            norm_values::Vector{Real} = Vector{Real}(undef, length(measured_elements))
            if lowercase(normalisation_values) == "CI_PO2016"
                for i in eachindex(measured_elements)
                    norm_values[i] = CI_CHONDRITE_PO2016[Symbol(measured_elements[i])][1]
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
        Ω::Diagonal{Float64x4,Vector{Float64x4}} = Diagonal(ω .^ 2)
        X̃::Matrix{Float64x4} = exp(-0.5log(Ω)) * X
        ỹ::Vector{Float64x4} = exp(-0.5log(Ω)) * y

        if cond(X̃) ≤ 1e7
            F = qr(X̃)
            Λ::Vector{Float64x4} = F \ ỹ
            VarΛX = Symmetric(inv(F.R) * transpose(inv(F.R)))
        else
            F = svd(X̃)
            Λ = F.V * inv(Diagonal(F.S)) * transpose(F.U) * ỹ
            VarΛX = F.V * inv(Diagonal(F.S .^2)) * F.Vt
        end
        rss::Vector{Float64x4} = Vector{Float64x4}(undef, 5)
        @simd for i ∈ eachindex(order)
            residuals::Vector{Float64x4} = (y .- (view(X, :, 1:i) * Λ[1:i]))
            rss[i] = transpose(residuals) * inv(Ω) * (residuals)
        end
        AIC::Vector{Float64x4} = Vector{Float64x4}(undef, 5)
        AIC = _akaike_information_criteria.(rss, 𝑁, order)

        for i in eachindex(Λ)
            Λ[i] = abs(Λ[i]) ≤ Base.rtoldefault(Float64x4) ? 0.0 : Λ[i]
        end
        mse = rss ./ (𝑁 .- (order .+ 1))
        Λ_SE::AbstractMatrix{Float64x4} = zeros(Float64x4, 5, 5)
        for i ∈ eachindex(order)
            Λ_SE[1:i, i] = sqrt.(diag(view(VarΛX, 1:i, 1:i) * (mse[i])))
        end
        tss::Float64x4 = transpose((y .- mean(y))) * inv(Ω) * (y .- mean(y))
        rmse::Vector{Float64x4} = sqrt.(mse)
        nrmse::Vector{Float64x4} = rmse ./ (maximum(y) - minimum(y))
        R²::Vector{Float64x4} = 1 .- (rss ./ (tss))
        R²ₒₚ::Vector{Float64x4} = _olkin_pratt.(R², 𝑁, order .+ 1)
        BIC::Vector{Float64x4} = Vector{Float64x4}(undef, 5)
        BIC = _bayesian_information_criteria.(rss, 𝑁, order)
        BICw =
            exp.(-0.5 .* (BIC .- minimum(BIC))) ./ sum(exp.(-0.5 .* (BIC .- minimum(BIC))))
        AIC = _akaike_information_criteria.(rss, 𝑁, order)
        AICw =
            exp.(-0.5 .* (AIC .- minimum(AIC))) ./ sum(exp.(-0.5 .* (AIC .- minimum(AIC))))
        return OrthogonalPolynomial(
            Float64.(Λ),
            UpperTriangular(Λ_SE),
            big.(β),
            big.(γ),
            big.(δ),
            big.(ϵ),
            Float64.(VarΛX),
            order,
            R²,
            R²ₒₚ,
            rmse,
            nrmse,
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
