#= Preamble
This file contains functions to load lanthanoid elemental data and calculate lambda fits.
=#

export REE_lambda

# Handler function for REE lambda fitting and anomaly calculation
function REE_lambda(
    data::DataFrame;
    measured_lanthanoids::Vector{String} = [
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
    ],
    uncertainties::Bool = false,
    uncertainty_type::AbstractString = "abs",
    uncertainty_level::Integer = 2,
    uncertainty_suffix::AbstractString = "_2SE",
    normalise::Bool = true,
    normalisation_values::AbstractString = "CI_PO2014",
    fit_ce::Bool = false,
    fit_eu::Bool = false,
    fit_gd::Bool = true,
    suppress_warnings::Bool = false,
)
    if uncertainties === false
        lanth_uncertainties = nothing
    else
        # measured_lanthanoids .* uncertainty_suffix
    end
    if nrow(data) > 1
        if !in("orthogonal_fit", names(data))
            insertcols!(
                data,
                :orthogonal_fit =>
                    Vector{GeochemistryTools.OrthogonalPolynomial}(undef, nrow(data)),
            )
        end
        X_REE = _X_REE()
        Threads.@threads for i ∈ 1:1:nrow(data)
            data.orthogonal_fit[i] = _lambdaREE(
                coalesce.(Vector(data[i, measured_lanthanoids]), NaN),
                lanth_uncertainties;
                lanth_measured = measured_lanthanoids,
                weight_type = uncertainty_type,
                normalise = normalise,
                normalisation_values = normalisation_values,
                fit_ce = fit_ce,
                fit_eu = fit_eu,
                fit_gd = fit_gd,
                operator_matrix = X_REE,
                suppress_warnings = suppress_warnings,
            )
        end
    end
    # _lambdaREE(Vector(select(data[i:])))
end

# Main function for calculating REE lambda
function _lambdaREE(
    lanth_values::AbstractVector{Float64},
    lanth_uncertainties::Union{Nothing,AbstractVector{T}} where {T<:Real} = nothing;
    lanth_measured::AbstractVector{String} = [
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
    ],
    weight_type::AbstractString = "abs",
    normalise::Bool = true,
    normalisation_values::AbstractString = "CI_PO2014",
    fit_ce::Bool = false,
    fit_eu::Bool = false,
    fit_gd::Bool = true,
    operator_matrix::Union{Nothing,Tuple} = nothing,
    suppress_warnings::Bool = false,
)
    fitting_indices::Vector{Int64} = findall(isfinite, lanth_values)
    if length(fitting_indices) ≤ 4
        if suppress_warnings != true
            @warn "Four or fewer values for sample, fit will be too unreliable."
        end
        return OrthogonalPolynomial(
            fill(nothing, length(fieldnames(OrthogonalPolynomial)))...,
        )
    else
        try
            if isnothing(operator_matrix)
                X, β, γ, δ, ϵ = _X_REE()
            else
                X, β, γ, δ, ϵ = operator_matrix
            end

            # if isnothing(lanth_uncertainties) != true
            #     fitting_indices =
            #         intersect(fitting_indices, findall(isfinite, lanth_uncertainties))
            #     lanth_uncertainties = lanth_uncertainties[fitting_indices]
            # end

            # remove declared anomalous elements
            if fit_gd == false
                deleteat!(
                    fitting_indices,
                    fitting_indices .== findfirst(occursin("Gd"), lanth_measured),
                )
            end
            if fit_eu == false
                deleteat!(
                    fitting_indices,
                    fitting_indices .== findfirst(occursin("Eu"), lanth_measured),
                )
            end
            if fit_ce == false
                deleteat!(
                    fitting_indices,
                    fitting_indices .== findfirst(occursin("Ce"), lanth_measured),
                )
            end

            # normalise data to reference values
            if normalise == true
                norm_values::MVector{length(fitting_indices),Float64} =
                    MVector{length(fitting_indices),Float64}(undef)
                if uppercase(normalisation_values) == "CI_PO2014"
                    lookup_dict = CI_CHONDRITE_PO2014
                elseif uppercase(normalisation_values) == "PAAS_P2012"
                    lookup_dict = PAAS_P2012
                end
                for i ∈ eachindex(norm_values)
                    norm_values[i] =
                        lookup_dict[Symbol(lanth_measured[fitting_indices][i])][1]
                end
                lanth_values = log.(lanth_values[fitting_indices] ./ norm_values)
            else
                lanth_values = log.(lanth_values[fitting_indices])
            end

            X = X[fitting_indices, :]

            order::SVector{5,Integer} = [0, 1, 2, 3, 4]
            if lanth_uncertainties === nothing
                ω::SVector{length(fitting_indices),Float64x4} =
                    fill(1.0, length(fitting_indices))
                # elseif occursin("rel", lowercase(weight_type)) === true
                #     ω = lanth_uncertainties[fitting_indices]
                # elseif occursin("abs", lowercase(weight_type)) == true
                #     ω =
                #         abs.(lanth_uncertainties[fitting_indices]) ./
                #         abs.(lanth_values[fitting_indices])
            else
                throw(
                    ArgumentError(
                        "Value of 'weight_type' is unrecognised. String should contain either 'rel' or 'abs'.",
                    ),
                )
            end
            𝑁::Integer = length(fitting_indices)
            Ω::Diagonal{Float64x4,Vector{Float64x4}} = Diagonal(ω .^ 2)
            X̃::Matrix{Float64x4} = exp(-0.5log(Ω)) * X
            ỹ::SVector{size(lanth_values)[1],Float64x4} = exp(-0.5log(Ω)) * lanth_values

            if cond(X̃) ≤ 1e7
                F = qr(X̃)
                Λ::MVector{5,Float64x4} = F \ ỹ
                VarΛX = Symmetric(inv(F.R) * transpose(inv(F.R)))
            else
                F = svd(X̃)
                Λ = F.V * inv(Diagonal(F.S)) * transpose(F.U) * ỹ
                VarΛX = F.V * inv(Diagonal(F.S .^ 2)) * F.Vt
            end
            rss::MVector{5,Float64x4} = MVector{5,Float64x4}(undef)
            @simd for i ∈ eachindex(order)
                residuals::MVector{size(lanth_values)[1],Float64x4} =
                    (lanth_values .- (view(X, :, 1:i) * Λ[1:i]))
                rss[i] = transpose(residuals) * inv(Ω) * (residuals)
            end

            for i ∈ eachindex(Λ)
                Λ[i] = abs(Λ[i]) ≤ Base.rtoldefault(Float64x4) ? 0.0 : Λ[i]
            end
            mse::SVector{5,Float64x4} = rss ./ (𝑁 .- (order .+ 1))
            Λ_SE::MMatrix{5,5,Float64x4} = zeros(Float64x4, 5, 5)
            for i ∈ eachindex(order)
                Λ_SE[1:i, i] = sqrt.(diag(view(VarΛX, 1:i, 1:i) * (mse[i])))
            end
            tss::Float64x4 =
                transpose((lanth_values .- mean(lanth_values))) *
                inv(Ω) *
                (lanth_values .- mean(lanth_values))
            rmse::SVector{5,Float64} = sqrt.(mse)
            nrmse::SVector{5,Float64} =
                rmse ./ (maximum(lanth_values) - minimum(lanth_values))
            R²::SVector{5,Float64} = 1 .- (rss ./ (tss))
            R²ₒₚ::SVector{5,Float64} = _olkin_pratt.(R², 𝑁, order .+ 1)
            AIC::SVector{5,Float64} = _akaike_information_criteria.(rss, 𝑁, order)
            BIC::SVector{5,Float64} = _bayesian_information_criteria.(rss, 𝑁, order)
            BICw::SVector{5,Float64} =
                exp.(-0.5 .* (BIC .- minimum(BIC))) ./
                sum(exp.(-0.5 .* (BIC .- minimum(BIC))))
            AICw::SVector{5,Float64} =
                exp.(-0.5 .* (AIC .- minimum(AIC))) ./
                sum(exp.(-0.5 .* (AIC .- minimum(AIC))))
            return OrthogonalPolynomial(
                Λ,
                UpperTriangular(Λ_SE),
                β,
                γ,
                δ,
                ϵ,
                VarΛX,
                order,
                R²,
                R²ₒₚ,
                rmse,
                nrmse,
                Float64.(rss),
                Float64.(mse),
                AIC,
                AICw,
                BIC,
                BICw,
                𝑁,
            )
        catch
            if suppress_warnings != true
                @warn "Error fitting model"
            end
            return OrthogonalPolynomial(
                fill(nothing, length(fieldnames(OrthogonalPolynomial)))...,
            )
        end
    end
end


function _X_REE(
    lanthanoids::Vector{String} = [
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
    ],
)
    lanth_radii = IONIC_RADIUS_CNEIGHT.(lanthanoids .* "_3+")
    𝑁::Integer = length(lanth_radii)
    x_sums::MVector{7,Float64x4} = Vector{Float64x4}(undef, 7)
    @simd for i ∈ eachindex(x_sums)
        x_sums[i] = sum(lanth_radii .^ i)
    end
    β::Float64x4 = _beta_orthogonal(𝑁, x_sums)
    γ::SVector{2,Float64x4} = _gamma_orthogonal(𝑁, x_sums)
    δ::SVector{3,Float64x4} = _delta_orthogonal(𝑁, x_sums)
    ϵ::SVector{4,Float64x4} = _epsilon_orthogonal(𝑁, x_sums)
    X::Matrix{Float64x4} = hcat(
        fill(1.0, 𝑁),
        (lanth_radii .- β),
        (lanth_radii .- γ[1]) .* (lanth_radii .- γ[2]),
        (lanth_radii .- δ[1]) .* (lanth_radii .- δ[2]) .* (lanth_radii .- δ[3]),
        (lanth_radii .- ϵ[1]) .* (lanth_radii .- ϵ[2]) .* (lanth_radii .- ϵ[3]) .*
        (lanth_radii .- ϵ[4]),
    )

    X = X[Not(findfirst(occursin("Pm"), lanthanoids)), :]

    return (X, β, γ, δ, ϵ)
end
