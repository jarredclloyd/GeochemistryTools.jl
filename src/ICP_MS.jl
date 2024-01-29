#=
This file contains tools for (Agilent) ICP-MS data
=#

export load_agilent, load_agilent2

"""
    load_agilent(args...; kwargs...)

Load and prepare data from CSV agilent CSV files.

Load all relevant files in the `host_directory` using glob, concatenates the data into one table, calculates ratio between `cps_column1` and `cps_column2`, the median-normalised ratio, and sorts the table by time.

# Arguments

  - `host_directory::AbstractString`: The folder path containing raw LA-ICP-MS CSV files.
  - `cps_column1::AbstractString`: A string representing the first desired column with counts per second data.
  - `cps_column2::AbstractString`: A string representing the second desired column with counts per second data.

## Argument Notes

  - `host_directory` will need to be wrapped by `raw` on Windows operating systems.
  - `cps_columnX` should be specified as `element` then `isotopic mass` (e.g., `Rb85`). It matches using regular
    expression and normalises the header names. If duplicate isotope values are present (due to mass shift checks, e.g.
    K39 -> 39, K39 -> 58) specify the desired mass by replacing `->` with `_` (e.g. `K39 -> 39` to `K39_39`).

# Keywords

  - `sample::AbstractString`: The sample name as a string: only needs to be a unique portion of the sample name. E.g. NIST610 will load all `NIST610*.csv` files, whereas NIST will load any file where the name contains `NIST*.csv`).
  - `date_time_format::AbstractString`: A string formatted for the appropriate DateTime in the CSV file. See `? Dates.DateFormat` for valid strings.
  - `first_row::Integer`: An `Integer` representing the first data row. Default value is `5`.
  - `gas_blank::Real`: The time (in seconds) where the gas blank ends. Default value is `27.5`.
  - `stable_time::Real`: The time (in seconds) where signal counts are stabilised. Default value is `32` (i.e. 32 seconds).
    This parameter is used to filter the rows when calculating the median counts per second ratio.
  - `signal_end::Real`: The time (in seconds) where the data should be truncated. Default value is `Inf` (will use entire
    signal beyond `stable_time`). It is recommended to set this to a value appropriate to your data (i.e. total signal
    time minus ~5 [e.g. 65]) to avoid some artefacts  that may occur near the end of an ablation.
  - `trim::Bool`: Whether to trim the data to only contain values between `stable_time` and `signal_end`. Default value is
    `false`. Setting to `true` it will retain only data for the time interval between `stable_time` and `signal_end`.

"""
function load_agilent(
    host_directory::AbstractString,
    cps_column1::AbstractString,
    cps_column2::AbstractString;
    sample::Union{Nothing,AbstractString} = nothing,
    date_time_format::AbstractString = "d/m/Y H:M:S",
    first_row::Integer = 5,
    gas_blank::Real = 27.5,
    stable_time::Real = 32,
    signal_end::Real = Inf,
    trim::Bool = false,
    aggregate_all::Bool = false,
    normalise::Bool = true,
    normalisation::AbstractString = "gmean"
)
    if aggregate_all === false && sample === nothing
        throw(
            ArgumentError(
                "Either `sample` needs to be specified or aggregate_all needs to be true",
            ),
        )
    end
    data = DataFrame()
    if aggregate_all === true
        files = glob("*.csv", host_directory)
    else
        files = glob(sample * "*.csv", host_directory)
    end
    for file in files
        df = CSV.read(
            file,
            DataFrame;
            header = false,
            limit = 3,
            silencewarnings = true,
            delim = ',',
        )
        analysis_name = chop(df[1, 1]; head = findlast("b\\", df[1, 1])[2], tail = 2)
        sample_name = chop(
            analysis_name;
            tail = length(analysis_name) - findlast("-", analysis_name)[1] + 1,
        )
        sample_name = rstrip(sample_name, ' ')
        analysis_time = chop(
            df[3, 1];
            head = findfirst(":", df[3, 1])[1] + 1,
            tail = length(df[3, 1]) -
                   findnext(" u", df[3, 1], findfirst(":", df[3, 1])[1])[1] + 1,
        )
        analysis_time = DateTime(analysis_time, date_time_format)
        if Dates.Year(analysis_time) < Dates.Year(2000)
            analysis_time = analysis_time + Dates.Year(2000)
        else
            analysis_time = analysis_time
        end
        df = CSV.read(
            file,
            DataFrame;
            header = 4,
            skipto = first_row,
            footerskip = 3,
            ignoreemptyrows = true,
            normalizenames = true,
            delim = ',',
        )
        df = select!(df, r"Time", r"" * cps_column1, r"" * cps_column2)
        rename!(df, ["time", cps_column1, cps_column2])
        insertcols!(df, 1, "sample" => sample_name)
        insertcols!(df, 2, "analysis_name" => analysis_name)
        insertcols!(df, 3, "analysis_time" => analysis_time)
        gas_blank_cps_column1 = geomean_zeros(df[0 .< df.time .< gas_blank, cps_column1])
        gas_blank_cps_column2 = geomean_zeros(df[0 .< df.time .< gas_blank, cps_column2])
        transform!(
            df,
            Cols(cps_column1, :time) =>
                ByRow(
                    (value, time) ->
                        time < mean([gas_blank, stable_time]) ? value : value - gas_blank_cps_column1,
                ) => cps_column1 * "_gbsub",
        )
        insertcols!(df, cps_column1 * "_σ" => sqrt.(abs.(df[!, cps_column1 * "_gbsub"])))
        transform!(
            df,
            Cols(cps_column2, :time) =>
                ByRow(
                    (value, time) ->
                        time < mean([gas_blank, stable_time]) ? value : value - gas_blank_cps_column2,
                ) => cps_column2 * "_gbsub",
        )
        insertcols!(df, cps_column2 * "_σ" => sqrt.(abs.(df[!, cps_column2 * "_gbsub"])))
        df.ratio = df[!, cps_column1 * "_gbsub"] ./ df[!, cps_column2 * "_gbsub"]
        insertcols!(
            df,
            "ratio_σ" =>
                sqrt.(
                    df[!, :ratio] .^ 2 .* (
                        (df[!, cps_column1 * "_σ"] ./ df[!, cps_column1 * "_gbsub"]) .^ 2 +
                        (df[!, cps_column2 * "_σ"] ./ df[!, cps_column2 * "_gbsub"]) .^ 2
                    )
                ),
        )
        if normalise == true
            if normalisation == "gmean"
                alg = geomean_zeros
            elseif normalisation == "median"
                alg = median
            elseif normalisation == "amean"
                alg = mean
            end
            df.ratio_norm =
            df.ratio ./ alg(df[stable_time .< df.time .< signal_end, :ratio])
            df.ratio_norm_σ = df.ratio_norm .* (df.ratio_σ ./ df.ratio)
        end
        if trim == true
            filter!(:time => x -> stable_time .< x .< signal_end, df)
        end
        append!(data, df)
    end
    sort!(data, :time)
    metadata!(data, "stable time", stable_time)
    metadata!(data, "signal end", signal_end)
    return data
end

"""
    load_agilent2(args...; kwargs...)

Load and prepare data from CSV agilent CSV files to assess downhole fractionation.

Load all relevant files in the `host_directory` using glob, concatenates the data into one table, calculates ratio between `cps_column1` and `cps_column2`, the median-normalised ratio, and sorts the table by time.

# Arguments

  - `host_directory::AbstractString`: The folder path containing raw LA-ICP-MS CSV files.
  - `sample::AbstractString`: The sample name as a string.
  - `cps_column1::AbstractString`: A string representing the first desired column with counts per second data.
  - `cps_column2::AbstractString`: A string representing the second desired column with counts per second data.

## Argument Notes

  - `host_directory` will need to be wrapped by `raw` on Windows operating systems.
  - `sample` only needs to be a unique portion of the sample name. E.g. NIST610 will just load `NIST610*.csv` files,
    whereas NIST will load any file where the name contains `NIST*.csv`).
  - `cps_columnX` should be specified as `element` then `isotopic mass` (e.g., `Rb85`). It matches using regular
    expression and normalises the header names. If duplicate isotope values are present (due to mass shift checks, e.g.
    K39 -> 39, K39 -> 58) specify the desired mass by replacing `->` with `_` (e.g. `K39 -> 39` to `K39_39`).

# Keywords

  - `first_row::Integer`: An `Integer` representing the first data row. Default value is `5`.
  - `gas_blank::Real`: The time (in seconds) where the gas blank ends. Default value is `27.5`.
  - `stable_time::Real`: The time (in seconds) where signal counts are stabilised. Default value is `32` (i.e. 32 seconds).
    This parameter is used to filter the rows when calculating the median counts per second ratio.
  - `signal_end::Real`: The time (in seconds) where the data should be truncated. Default value is `Inf` (will use entire
    signal beyond `stable_time`). It is recommended to set this to a value appropriate to your data (i.e. total signal
    time minus ~5 [e.g. 65]) to avoid some artefacts  that may occur near the end of an ablation.
  - `trim::Bool`: Whether to trim the data to only contain values between `first_row` and `signal_end`. Default value is
    `false`. Setting to `true` it will delete all data from the first record past `signal_end` to the last record.
"""
function load_agilent2(
    host_directory::AbstractString,;
    sample::Union{Nothing,AbstractString} = nothing,
    date_time_format::AbstractString = "d/m/Y H:M:S",
    first_row::Integer = 5,
    gas_blank::Real = 27.5,
    stable_time::Real = 32,
    signal_end::Real = Inf,
    trim::Bool = false,
    aggregate_all::Bool = false,
)
    if aggregate_all === false && sample === nothing
        throw(
            ArgumentError(
                "Either `sample` needs to be specified or aggregate_all needs to be true",
            ),
        )
    end
    data = DataFrame()
    if aggregate_all === true
        files = glob("*.csv", host_directory)
    else
        files = glob(sample * "*.csv", host_directory)
    end
    for file in files
        df = CSV.read(
            file,
            DataFrame;
            header = false,
            limit = 3,
            silencewarnings = true,
            delim = ',',
        )
        analysis_name = chop(df[1, 1]; head = findlast("b\\", df[1, 1])[2], tail = 2)
        sample_name = chop(
            analysis_name;
            tail = length(analysis_name) - findlast("-", analysis_name)[1] + 1,
        )
        sample_name = rstrip(sample_name, ' ')
        analysis_time = chop(
            df[3, 1];
            head = findfirst(":", df[3, 1])[1] + 1,
            tail = length(df[3, 1]) -
                   findnext(" u", df[3, 1], findfirst(":", df[3, 1])[1])[1] + 1,
        )
        analysis_time = DateTime(analysis_time, date_time_format)
        if Dates.Year(analysis_time) < Dates.Year(2000)
            analysis_time = analysis_time + Dates.Year(2000)
        else
            analysis_time = analysis_time
        end
        df = CSV.read(
            file,
            DataFrame;
            header = 4,
            skipto = first_row,
            footerskip = 3,
            ignoreemptyrows = true,
            normalizenames = true,
            delim = ',',
        )
        rename!(df, "Time_Sec_" => "time")
        insertcols!(df, 1, "sample" => sample_name)
        insertcols!(df, 2, "analysis_name" => analysis_name)
        insertcols!(df, 3, "analysis_time" => analysis_time)
        append!(data, df)
    end
    metadata!(data, "stable time", stable_time)
    metadata!(data, "signal end", signal_end)
    return data
end


function automatic_laser_times(
    time,
    signal;
    bandwidth = ceil(Integer, sqrt(length(signal))/2),
)
    medians = Vector{Float64}(undef, cld(length(signal), bandwidth))
    for i in eachindex(medians)
        medians[i] = median(signal[round(Int, max((i - 1) * bandwidth, firstindex(signal))):round(Int, min(i * bandwidth, lastindex(signal)))])
    end
    pvalue_KS = pvalue(ApproximateTwoSampleKSTest(medians[begin: round(Int, end / 2)], medians[round(Int, end / 2):end]))
    pvalue_JB = pvalue(JarqueBeraTest(signal))
    if pvalue_KS > 0.05 || pvalue_JB > 0.05
        println("no signal detected")
    else
        z = Array{Float64}(undef, length(signal), 5)
        for i in eachindex(signal)
            z[i, 1] = if i < bandwidth
                geomean_zeros(signal[begin:i])
            else
                geomean_zeros(signal[(i - bandwidth + 1):i])
            end
            z[i, 2] = if i < bandwidth
                geovar_zeros(signal[begin:i])
            else
                geovar_zeros(signal[(i - bandwidth + 1):i])
            end
            z[i, 3] = if i < bandwidth
                -(z[i, 1], z[begin, 1]) / -(time[i], time[begin])
            else
                -(z[i, 1], z[i - bandwidth + 1, 1]) / -(time[i], time[i - bandwidth + 1])
            end
            z[i, 4] = i < 3 ? -Inf : -(z[i, 3], z[i - 1, 3])
            z[i, 5] = i < 3 ? Inf : pvalue(OneSampleTTest(z[2:i, 2], mean(z[2:(i - 1), 2])))
        end
        z[1, 2:3] .= -Inf
        laser_start = time[findmin(z[:, 5])[2]]
        aerosol_arrival = time[findmax(z[:, 4])[2]]
        return (laser_start, aerosol_arrival, z) # remove z when algorithm finalised
    end
end
