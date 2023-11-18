export fit_eivlr, affine_prediction, affine_standarderror, affine_confidenceband, affine_predictionband
export LinearRegression, ErrorsInVariablesRegression
abstract type LinearRegression <: Any end
abstract type ErrorsInVariablesRegression <: LinearRegression end

Base.getindex(LR::LinearRegression, i::Integer) = getfield(LR, i)
Base.keys(LR::LinearRegression) = Base.keys(LR)

struct York <: ErrorsInVariablesRegression
    beta0::Real
    beta0_se::Real
    beta1::Real
    beta1_se::Real
    x_intercept::Real
    x_intercept_se::Real
    reduced_chi_squared::Real
    p_value::Real
    covariance_beta::Real
    n_observations::Integer
    x_bar::Real
    y_bar::Real
end

function Base.show(io::IOContext, F::York)
    println(
        io,
        "βₒ: $(round(F.beta0, sigdigits = 5)) ± $(round(F.beta0_se, sigdigits = 5))",
    )
    println(
        io,
        "β₁: $(round(F.beta1, sigdigits = 5)) ± $(round(F.beta1_se, sigdigits = 5))",
    )
    println(
        io,
        "x_intercept: $(round(F.x_intercept, sigdigits = 5)) ± $(round(F.x_intercept_se, sigdigits = 5))",
    )
    println(
        io,
        "χ²ᵣ (MSWD): $(round(F.reduced_chi_squared, sigdigits = 3)); p-value: $(round(F.p_value, sigdigits = 5))",
    )
    println(io, "no. observations: $(F.n_observations)")
end

struct MahonNonFixed <: ErrorsInVariablesRegression
    beta0::Real
    beta0_se::Real
    beta1::Real
    beta1_se::Real
    x_intercept::Real
    x_intercept_se::Real
    reduced_chi_squared::Real
    p_value::Real
    covariance_beta::Real
    n_observations::Integer
end

function Base.show(io::IOContext, F::MahonNonFixed)
    println(
        io,
        "βₒ: $(round(F.beta0, sigdigits = 5)) ± $(round(F.beta0_se, sigdigits = 5))",
    )
    println(
        io,
        "β₁: $(round(F.beta1, sigdigits = 5)) ± $(round(F.beta1_se, sigdigits = 5))",
    )
    println(
        io,
        "x_intercept: $(round(F.x_intercept, sigdigits = 5)) ± $(round(F.x_intercept_se, sigdigits = 5))",
    )
    println(
        io,
        "χ²ᵣ (MSWD): $(round(F.reduced_chi_squared, sigdigits = 3)); p-value: $(round(F.p_value, sigdigits = 5))",
    )
    println(io, "no. observations: $(F.n_observations)")
end

struct MahonFixed <: ErrorsInVariablesRegression
    beta0::Real
    beta0_se::Real
    beta1::Real
    beta1_se::Real
    x_intercept::Real
    x_intercept_se::Real
    reduced_chi_squared::Real
    p_value::Real
    covariance_beta::Real
    n_observations::Integer
    fixed_point::Tuple{Real,Real,Real,Real}
end

function Base.show(io::IOContext, F::MahonFixed)
    println(
        io,
        "βₒ: $(round(F.beta0, sigdigits = 5)) ± $(round(F.beta0_se, sigdigits = 5))",
    )
    println(
        io,
        "β₁: $(round(F.beta1, sigdigits = 5)) ± $(round(F.beta1_se, sigdigits = 5))",
    )
    println(
        io,
        "x_intercept: $(round(F.x_intercept, sigdigits = 5)) ± $(round(F.x_intercept_se, sigdigits = 5))",
    )
    println(
        io,
        "χ²ᵣ (MSWD): $(round(F.reduced_chi_squared, sigdigits = 3)); p-value: $(round(F.p_value, sigdigits = 5))",
    )
    println(io, "no. observations: $(F.n_observations)")
    println(io, "fixed point: $(F.fixed_point)")
end

struct GeneralisedLeastSquares <: LinearRegression
    beta::AbstractVector
    beta_se::AbstractVector
    variance_covariance::AbstractArray
    r_squared::AbstractFloat
    rmse::AbstractFloat
    n_observations::Integer
end

function Base.show(io::IOContext, F::GeneralisedLeastSquares)
    for i in eachindex(F.btea)
        println(
            io,
            "β$(i - 1): $(round.(F.beta[i]; digits = 5)) ± $(round.(F.beta_se[i]; digits = 5))",
        )
    end
    println(io, "R²: $(round(F.r_squared; digits = 4))")
    println(io, "RMSE = $(round(F.rmse; digits=4))")
    println(io, "no. observations: $(F.n_observations)")
end

"""
    fit_eivlr(df, algorithm="mahon"; kwargs..)

Compute line of best fit using an errors-in-variables linear regression algorithm.

Input df as a DataFrame of 4 of 5 columns wide with column order (X, sX, Y, sY, [ρXY]).
Algorithms available are '"mahon"' and '"york"'.

# Keywords

  - `se_level_in::Integer`: Standard error level of input data. Provide as a positive integer.
  - `se_type::AbstractString`: Standard error type as a string of value `"abs"` OR `"rel"`. Values equal to
    `'a'`, `"absolute"`, `'r'`, and `"relative"` will also work. Case insensitive.
  - `initial::Any`: A value for the y-intercept. Can be input as a string key from an appropriate dictionary, as a single
    numeric value, or as a vector of the initial and its standard error (same `se_level_in` as input data). E.g. initial =
    "MDCInv", initial = 0.72, OR initial = [0.72, 0.01].

      + Dictionaries available are `dict_sr87_sr86i`
      + For a full list of available keys in any dictionary type `keys(<dict_name>)`

# Example

```julia-repl
julia> fit_eivlr(df, "mahon"; se_level_in = 2, se_type = "abs", initial = "MDCInv")

```

# References

York D. et al. 2004 "Unified equations for the slope, intercept, and standard errors of the best straight line", *American
Journal of Physics*, 72(3), doi:https://doi.org/10.1119/1.1632486.

Mahon, KI (1996) 'The New “York” Regression: Application of an Improved Statistical Method to Geochemistry',
*International Geology Review*, 38(4):293–303. https://doi.org/10.1080/00206819709465336

Stephan, T & Trappitsch, R (2023) 'Reliable uncertainties: Error correlation, rotated error bars, and linear
regressions ∈ three-isotope plots and beyond', *International Journal of Mass Spectrometry*, 491:117053.
https://doi.org/10.1016/j.ijms.2023]].117053
"""
function fit_eivlr(
    df::AbstractDataFrame,
    algorithm::AbstractString = "mahon";
    se_level_in::Integer = 2,
    se_type::AbstractString = "abs",
    initial::Any = nothing,
)
    dfCols::Integer = ncol(df)
    dfRows::Integer = nrow(df)
    if initial !== nothing
        if isa(initial, AbstractString) == true && haskey(dict_sr87_sr86i, initial) == true
            if algorithm .== "mahon"
                initial = [
                    0,
                    0,
                    deepcopy(dict_sr87_sr86i[initial])[1],
                    deepcopy(dict_sr87_sr86i[initial])[2] * se_level_in,
                ]
            else
                initial = [
                    1e-16,
                    1e-16,
                    deepcopy(dict_sr87_sr86i[initial])[1],
                    deepcopy(dict_sr87_sr86i[initial])[2] * se_level_in,
                ]
            end
        elseif length(initial) == 1 && isa(initial, AbstractFloat)
            if algorithm .== "mahon"
                initial = [0, 0, initial, initial * 0.1]
            else
                initial = [1e-16, 1e-16, initial, initial * 0.1]
            end
        elseif length(initial) == 2 &&
               isa(initial[1], AbstractFloat) &&
               isa(initial[2], AbstractFloat)
            if algorithm .== "mahon"
                initial = [0, 0, initial[1], initial[2]]
            else
                initial = [1e-16, 1e-16, initial[1], initial[2]]
            end
        elseif length(initial) >= 3
            throw(
                ArgumentError(
                    "Only the initial ratio value and its uncertainty is required.",
                ),
            )
        else
            throw(
                ArgumentError(
                    """
initial variable is malformed.
Use either a string that matches a fixed value or reference material name:
    e.g. "MDC", "MDCInv", or "SolarSystem"
Or input a single numeric value or two-length numeric vector:
    e.g. 0.72 or [0.7200, 0.0010]

For a full list of available reference material keys in a relevant dictionary: keys(<dict_name>)

Dictionaries currently available are `dict_sr87_sr86i`
""",
                ),
            )
        end
        if dfCols == 4
            push!(df, initial)
        elseif dfCols == 5
            push!(initial, 0)
            push!(df, initial)
        end
    end

    if lowercase.(algorithm) == "mahon"
        alg = _eivlr_mahon
    elseif lowercase.(algorithm) == "york"
        alg = _eivlr_york
    end
    if occursin('a', lowercase.(se_type)) == true ||
       occursin("abs", lowercase.(se_type)) == true ||
       occursin("absolute", lowercase.(se_type)) == true
        if dfCols == 5
            eivlr = alg(
                df[!, 1],
                df[!, 2] ./ se_level_in,
                df[!, 3],
                df[!, 4] ./ se_level_in,
                df[!, 5],
            )
        elseif dfCols == 4
            eivlr = alg(df[!, 1], df[!, 2] ./ se_level_in, df[!, 3], df[!, 4] ./ se_level_in)
        else
            throw(
                ArgumentError("Column width is not equal to 4 or 5. Some data is missing."),
            )
        end
    elseif occursin('r', lowercase.(se_type)) == true ||
           occursin("rel", lowercase.(se_type)) == true ||
           occursin("relative", lowercase.(se_type)) == true
        if dfCols == 5
            eivlr = alg(
                df[!, 1],
                (df[!, 2] .* df[!, 1]) ./ se_level_in,
                df[!, 3],
                (df[!, 4] .* df[!, 3]) ./ se_level_in,
                df[!, 5],
            )
        elseif dfCols == 4
            eivlr = alg(
                df[!, 1],
                (df[!, 2] .* df[!, 1]) ./ se_level_in,
                df[!, 3],
                (df[!, 4] .* df[!, 3]) ./ se_level_in,
            )
        else
            throw(
                ArgumentError("Column width is not equal to 4 or 5. Some data is missing."),
            )
        end
    end
    if initial !== nothing
        popat!(df, dfRows + 1)
    end
    return eivlr
end

function affine_prediction(x::AbstractVector, fit::ErrorsInVariablesRegression)
    return fit.beta0 .+ x .* fit.beta1
end

function affine_standarderror(
    x::AbstractVector,
    fit::ErrorsInVariablesRegression;
    se_level::Integer = 1,
)
    return vec(sqrt.(abs.(fit.reduced_chi_squared .* (x * fit.covariance_beta)))) .* se_level
end
function affine_confidenceband(
    x::AbstractVector,
    fit::ErrorsInVariablesRegression;
    confidence_level::AbstractFloat = 0.95,
)
    ν = fit.n_observations - 2
    tvalue = cquantile(TDist(ν), (1 - confidence_level) / 2)
    return vec(sqrt.(abs.(fit.reduced_chi_squared .* (x * fit.covariance_beta)))) .* tvalue
end

function affine_predictionband(
    x::AbstractVector,
    fit::ErrorsInVariablesRegression;
    confidence_level::AbstractFloat = 0.95,
)
    ν = fit.n_observations - 2
    tvalue = cquantile(TDist(ν), (1 - confidence_level) / 2)
    return vec(sqrt.(abs.(fit.reduced_chi_squared .* (1 .+ x .* fit.covariance_beta)))) .*
           tvalue
end
