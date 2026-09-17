# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

export UNCERT_REG, uncertainty_source, uncertainty_report, uncertainty_contributions, reset!

export @uncertainty_source

# Types

# Uncertainty Report
struct UncertaintyReport
    result_name::String
    result::Measurement
    names::Vector{String}
    ids::Vector{UInt64}
    values::Vector{Float64}
    uncertainties::Vector{Float64}
    sensitivity_coeffs::Vector{Float64}
    contributions::Vector{Float64}
    fractional_contributions::Vector{Float64}
end

Base.getindex(uncert_rep::UncertaintyReport, key::Symbol) = getfield(uncert_rep, key)
function Base.getindex(uncert_rep::UncertaintyReport, key::AbstractString)
    return getfield(uncert_rep, Symbol(key))
end
function Base.haskey(uncert_rep::UncertaintyReport, key::Symbol)
    return key ∈ fieldnames(typeof(uncert_rep))
end
Base.keys(uncert_rep::UncertaintyReport) = fieldnames(typeof(uncert_rep))
Base.values(uncert_rep::UncertaintyReport) = getfield.(Ref(uncert_rep), keys(uncert_rep))
Base.pairs(uncert_rep::UncertaintyReport) = zip(keys(uncert_rep), values(uncert_rep))

const UNCERT_REG = Dict{Tuple{Float64,Float64,UInt},String}()

reset!() = empty!(UNCERT_REG)

function with_registry(f)
    reset!()
    try
        return f()
    finally
        reset!()
    end
end

macro uncertainty_source(ex)
    if ex isa Expr && ex.head == :(=)
        # case: @source x = measurement(...)
        var = ex.args[1]
        val = ex.args[2]
        return quote
            $(esc(var)) = $(esc(val))
            register_source!($(esc(var)), $(QuoteNode(var)))
        end
    else
        # case: @source x
        return :(register_source!($(esc(ex)), $(QuoteNode(ex))))
    end
end

# Helpers
function uncertainty_source(name::Union{String,Symbol}, meas::Measurement)
    return register_source!(meas, name)
end

function uncertainty_source(name::Union{String,Symbol}, val::Real, unc::Real)
    return register_source!(measurement(val, unc), name)
end

function register_source!(meas::Measurement, name::Union{String,Symbol})
    comps = Measurements.uncertainty_components(meas)
    length(comps) == 1 || error("Expected one independent uncertainty_source")

    push!(UNCERT_REG, only(keys(comps)) => String(name))

    return meas
end

# core function
function uncertainty_report(meas::Measurement; result_name::String = "Result")
    meas_var = Measurements.uncertainty(meas)^2
    comps = Measurements.uncertainty_components(meas)
    isempty(comps) && error("No uncertainty components")

    ks = collect(keys(comps))
    𝑛 = length(ks)

    𝑥 = Vector{Float64}(undef, 𝑛)
    𝑢 = Vector{Float64}(undef, 𝑛)
    𝑢ᵢ = Vector{Float64}(undef, 𝑛)
    𝑐ᵢ = Vector{Float64}(undef, 𝑛)
    names = Vector{String}(undef, 𝑛)
    ids = Vector{UInt64}(undef, 𝑛)
    frac_𝑢ᵢ = Vector{Float64}(undef, 𝑛)

    for (i, k) ∈ enumerate(ks)
        𝑥[i] = k[1]
        𝑢[i] = k[2]
        ids[i] = k[3]
        𝑐ᵢ[i] = meas.der[k]
        𝑢ᵢ[i] = comps[k]
        frac_𝑢ᵢ[i] = 𝑢ᵢ[i]^2 / meas_var
        names[i] = get(UNCERT_REG, k, string(k))
    end

    return UncertaintyReport(result_name, meas, names, ids, 𝑥, 𝑢, 𝑐ᵢ, 𝑢ᵢ, frac_𝑢ᵢ)
end

# Table builder for export
function uncertainty_contributions(rep::UncertaintyReport)

    u2y = Measurements.uncertainty(rep.result)^2
    u2y == 0 && return DataFrame()

    rows = NamedTuple[]

    for i ∈ eachindex(rep.names)
        var_i = rep.contributions[i]^2
        frac  = var_i / u2y

        push!(
            rows,
            (
                result = rep.result_name,
                input = rep.names[i],
                x_i = rep.values[i],
                u_xi = rep.uncertainties[i],
                u_i = rep.contributions[i],
                fraction = frac
            )
        )
    end

    df = DataFrame(rows)
    sort!(df, :fraction; rev = true)

    return df
end
