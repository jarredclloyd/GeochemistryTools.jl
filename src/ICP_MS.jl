#=
This file contains tools for (Agilent) ICP-MS data
=#

export load_agilent, load_agilent2, automatic_laser_times

"""
    load_agilent(args...; kwargs...)

Load and prepare data from CSV agilent CSV files.

Load all relevant files in the `host_directory` using glob, concatenates the data into one table, calculates ratio between `cps_column1` and `cps_column2`, the centred ratio, and sorts the table by time.

# Arguments

  - `host_directory::AbstractString`: The folder path containing raw LA-ICP-MS CSV files.
  - `cps_column1::AbstractString`: A string representing the first desired column with counts per second data.
  - `cps_column2::AbstractString`: A string representing the second desired column with counts per second data.

## Argument Notes

  - `host_directory` will need to be wrapped by `raw` on Windows operating systems.
  - `cps_columnX` should be specified as `element` then `isotopic mass` (e.g., `Rb85`). It matches using regular
    expression and centres the header names. If duplicate isotope values are present (due to mass shift checks, e.g.
    K39 -> 39, K39 -> 58) specify the desired mass by replacing `->` with `_` (e.g. `K39 -> 39` to `K39_39`).

# Keywords

  - `sample::AbstractString`: The sample name as a string: only needs to be a unique portion of the sample name. E.g. NIST610 will load all `NIST610*.csv` files, whereas NIST will load any file where the name contains `NIST*.csv`).
  - `date_time_format::AbstractString`: A string formatted for the appropriate DateTime in the CSV file. See `? Dates.DateFormat` for valid strings.
  - `first_row::Integer`: An `Integer` representing the first data row. Default value is `5`.
  - `automatic_times::Bool`: `[true | false]` Should laser, gas blank, and signal times be automatically determined. Default is `true`.
  - `gas_blank::Real`: The time (in seconds) where the gas blank ends. Default value is `27.5`. Is only used if `automatic_times == false`
  - `stable_time::Real`: The time (in seconds) where signal counts are stabilised. Default value is `32` (i.e. 32 seconds).
    This parameter is used to filter the rows when calculating the central tendency for counts per second ratio. Is only used if `automatic_times == false`
  - `signal_end::Real`: The time (in seconds) where the data should be truncated. Default value is `Inf` (will use entire
    signal beyond `stable_time`). It is recommended to set this to a value appropriate to your data (i.e. total signal
    time minus ~5 [e.g. 65]) to avoid some artefacts  that may occur near the end of an ablation. Is only used if `automatic_times == false`
  - `trim::Bool`: `[true | false]` Whether to trim the data to only contain values between `stable_time` and `signal_end`. Default value is
    `false`. Setting to `true` it will retain only data for the time interval between `stable_time` and `signal_end`. `automatic_times == true` will set these values.
  - `aggregate_files:Bool`: `[true | false]` A boolean flag for concatenating all relevant files in a directory. Default is `false`
  - `centre::Bool`: `[true | false]` A boolean flag to centre data around the `central_tendency` (i.e. shift to mean of zero). Default is `true`.
  - `central_tendency::AbstractString`: `["amean" | "gmean" | "median"]` The algorithm to use for the measure of central tendency (e.g. arithmetic mean, geometric mean, median)

"""
function load_agilent(
    host_directory::AbstractString,
    cps_column1::AbstractString,
    cps_column2::AbstractString;
    sample::Union{Nothing,AbstractString} = nothing,
    date_time_format::Union{AbstractString, DateFormat} = DateFormat("d/m/Y H:M:S"),
    first_row::Integer = 5,
    automatic_times::Bool = true,
    gas_blank::Real = 27.5,
    stable_time::Real = 32,
    signal_end::Real = Inf,
    trim::Bool = false,
    aggregate_files::Bool = false,
    centre::Bool = true,
    central_tendency::AbstractString = "gmean",
    spot_size::Union{Missing, AbstractString, Integer} = missing,
)
    if aggregate_files === false && sample === nothing
        throw(
            ArgumentError(
                "Either `sample` needs to be specified or aggregate_files needs to be true",
            ),
        )
    end
    if typeof(date_time_format) <: AbstractString
        date_time_format = Dates.DateFormat(date_time_format)
    end
    data = DataFrame()
    if aggregate_files === true
        files = glob("*.csv", host_directory)
    else
        files = glob(sample * "*.csv", host_directory)
    end
    for file in files
        head_info = split(readuntil(file, "Time "), "\n")
        analysis_name = chop(head_info[1]; head = findlast("\\", head_info[1])[1], tail = 3)
        sample_name = chop(
            analysis_name;
            tail = length(analysis_name) - findlast("-", analysis_name)[1] + 1,
        )
        sample_name = rstrip(sample_name, ' ')
        if ismissing(spot_size) !== true
            if typeof(spot_size) <: AbstractString
                spot_size = tryparse(Int, file[(findlast("_", file)[1] + 1):(findnext(
                    " ",
                    file,
                    findlast("_", file)[1],
                )[1] - 1)])
            elseif typeof(spot_size) <: Number
            else
                error("malformed spot_size variable")
            end
        end
        analysis_time = chop(
            head_info[3][findfirst(":",  head_info[3])[1] + 2 : findlast(":",  head_info[3])[1] + 2],
        )
        analysis_time = DateTime(analysis_time, date_time_format)
        if Dates.Year(analysis_time) < Dates.Year(2000)
            analysis_time = analysis_time + Dates.Year(2000)
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
        if automatic_times == true
            transform!(df, AsTable(Not(1)) => ByRow(sum) => :total_signal)
            auto_times = automatic_laser_times(df[!, 1], df[!, :total_signal])
            if !isnothing(auto_times) == true
                gas_blank_ind = auto_times[1][1]
                laser_time = auto_times[2][2]
                stable_time = auto_times[4][2]
                signal_end = auto_times[5][2]
            end
        else
            gas_blank_ind = findlast(≤(gas_blank), df[!, 1])
            laser_time = middle(gas_blank, stable_time)
        end
        if automatic_times == true && isnothing(auto_times) == true
            df = nothing
        else
            df = select!(df, r"Time", r"" * cps_column1, r"" * cps_column2)
            rename!(df, ["signal_time", cps_column1, cps_column2])
            insertcols!(df, 1, "sample" => sample_name)
            insertcols!(df, 2, "analysis_name" => analysis_name)
            insertcols!(df, 3, "spot_size" => spot_size)
            insertcols!(df, 4, "analysis_time" => analysis_time)
            insertcols!(df, 6, "beam_time" => df.signal_time .- laser_time)
            gas_blank_cps_column1 = geomean_zeros(df[begin:gas_blank_ind, cps_column1])
            gas_blank_cps_column2 = geomean_zeros(df[begin:gas_blank_ind, cps_column2])
            transform!(
                df,
                Cols(cps_column1, :signal_time) =>
                    ByRow(
                        (value, time) ->
                            time ≤ laser_time ? value : value - gas_blank_cps_column1,
                    ) => cps_column1 * "_gbsub",
            )
            insertcols!(
                df,
                cps_column1 * "_σ" => sqrt.(abs.(df[!, cps_column1 * "_gbsub"])),
            )
            transform!(
                df,
                Cols(cps_column2, :signal_time) =>
                    ByRow(
                        (value, time) ->
                            time ≤ laser_time ? value : value - gas_blank_cps_column2,
                    ) => cps_column2 * "_gbsub",
            )
            insertcols!(
                df,
                cps_column2 * "_σ" => sqrt.(abs.(df[!, cps_column2 * "_gbsub"])),
            )
            transform!(
                df,
                Cols(cps_column1 * "_gbsub", cps_column2 * "_gbsub", :signal_time) =>
                    ByRow(
                        (cps1, cps2, time) ->
                            time ≤ laser_time && iszero(cps2) == true ? 0.0 : cps1 / cps2,
                    ) => :ratio,
            )
            insertcols!(
                df,
                "ratio_σ" =>
                    sqrt.(
                        df[!, :ratio] .^ 2 .* (
                            (df[!, cps_column1 * "_σ"] ./ df[!, cps_column1 * "_gbsub"]) .^
                            2 +
                            (df[!, cps_column2 * "_σ"] ./ df[!, cps_column2 * "_gbsub"]) .^
                            2
                        )
                    ),
            )
            if centre == true
                if central_tendency == "gmean"
                    alg = geomean_zeros
                elseif central_tendency == "median"
                    alg = median
                elseif central_tendency == "amean"
                    alg = mean
                end
                centre_value = alg(df[stable_time .≤ df.signal_time .≤ signal_end, :ratio])
                df.ratio_centred =
                    df.ratio .- centre_value
                df.ratio_centred_σ = df.ratio_centred .* (df.ratio_σ ./ df.ratio)
            end
            if trim == true
                filter!(:signal_time => x -> stable_time .≤ x .≤ signal_end, df)
            end
            append!(data, df)
        end
    end
    if !isempty(data) == true
        sort!(data, :signal_time)
        return data
    end
end

"""
    load_agilent2(host_directory::AbstractString; kwargs...)

Load and prepare data from CSV agilent CSV files to assess downhole fractionation.

Load all relevant files in the `host_directory` using glob, concatenates the data into one
table, the median-centred ratio, and sorts the table by time.

# Arguments

  - `host_directory::AbstractString`: The folder path containing raw LA-ICP-MS CSV files.

# Keywords

  - `sample::AbstractString`: The sample name as a string.
  - `header_row::Integer`: An `Integer` representing the header row. Default value is `4`.
  - `first_row::Integer`: An `Integer` representing the first data row. Default value is `5`.
  - `footer_skip::Integer`: An `Integer` representing the number of trailing rows to trim. Default value is `3`.
  - `aggregate_files::Bool`: Default `false`. Will concatenate all CSV files in host folder to one data frame.

## Argument Notes

  - `host_directory` will need to be wrapped by `raw` on Windows operating systems if not using `/` or `\\\\`.

  - `date_time_format::AbstractString`: Default "d/m/Y H:M:S". A string formatted for the appropriate DateTime in the CSV file.

      + See `? Dates.DateFormat` for valid strings.
  - If `aggregate_files = true`, `sample` only needs to be a unique portion of the sample
    name. E.g. "NIST610" will load only `NIST610*.csv` files, whereas "NIST" will load any
    file where the name contains `NIST*.csv`).
  - If `aggregate_files = false`, `sample` needs to be the full sample name, e.g. "NIST610 - 1".
  - Either `sample` or `aggregate_files=true` need to be specified for the function to run.
"""
function load_agilent2(
    host_directory::AbstractString,
    ;
    sample::Union{Nothing,AbstractString} = nothing,
    date_time_format::AbstractString = "d/m/Y H:M:S",
    header_row::Integer = 4,
    first_row::Integer = 5,
    footer_skip::Integer = 3,
    aggregate_files::Bool = false,
)
    if aggregate_files === false && sample === nothing
        throw(
            ArgumentError(
                "Either `sample` needs to be specified or aggregate_files needs to be true",
            ),
        )
    end
    data = DataFrame()
    if aggregate_files === true && sample === nothing
        files = glob("*.csv", host_directory)
    elseif aggregate_files === false && sample !== nothing
        files = glob(sample * ".csv", host_directory)
    elseif aggregate_files === true && sample !== nothing
        files = glob(sample * "*.csv", host_directory)
    end
    for file in files
        head_info = split(readuntil(file, "Time "), "\n")
        analysis_name = chop(head_info[1]; head = findlast("\\", head_info[1])[1], tail = 3)
        sample_name = chop(
            analysis_name;
            tail = length(analysis_name) - findlast("-", analysis_name)[1] + 1,
        )
        sample_name = rstrip(sample_name, ' ')
        analysis_time = chop(
            head_info[3][(findfirst(":", head_info[3])[1] + 2):(findlast(
                ":",
                head_info[3],
            )[1] + 2)],
        )
        analysis_time = DateTime(analysis_time, date_time_format)
        if Dates.Year(analysis_time) < Dates.Year(2000)
            analysis_time = analysis_time + Dates.Year(2000)
        end
        df = CSV.read(
            file,
            DataFrame;
            header = header_row,
            skipto = first_row,
            footerskip = footer_skip,
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
    return data
end

"""
    automatic_laser_times(time::AbstractVector{<:Real}, signal::AbstractVector{<:Real}; kwargs...)

Attempt to automatically determine gas blank period, laser start, aerosol arrival time,
and signal period.

Is most accurate when using the total signal sum, but can work on single mass values.

# Extended Description

To determine if any signal present, the algorithm uses both a Jarque Bera Test and median
filter over the bandwidth to perform a two sample Kolmogorov-Smirnov test.
If either test returns a pvalue > 0.05, no signal is detected and the function stops.

If signal is determined to be present, smooths the signal using a geometric mean of the
signal over the bandwidth, and computes the geometric variance of the signal over the
bandwidth. Then uses a one sample T-test of the variances to determine the laser start
(minimum pvalue).

Then uses the upper and lower 0.05 quantiles of the smoothed signal from after the laser
start to find the aerosol arrival time and signal times.

If signal end time and signal start time are within 5 seconds of each other, algorithm will
attempt to re-estimate signal times determination if slope > 0. May produce inaccurate
signal times in this case if there is a small amount of signal (i.e. < 5 seconds).

Returns a tuple of tuples, (ind, time) for gas blank (end), laser start, aerosol arrival,
signal start & signal end. Likely to be updated to a `struct` in future.

# Keywords

  - `bandwidth::Integer`: Controls width of smoothing. Default is `cld(√(length(signal)), 2)`
    (this value generally works very well). Values too large will cause algorithm to
    fail to detect signal. Values too small may not accurately determine laser start.
  - `verbose::Bool`: Default is `false`. If set to `true` will print estimated times and
    indexes to REPL.

# Example

```julia-repl
julia> automatic_laser_times(G00.time, G00.signal_total)
((63, 28.6888), (65, 29.5973), (70, 31.8685), (71, 32.3228), (153, 69.5703))
```
"""
function automatic_laser_times(
    time::AbstractVector{<:Real},
    signal::AbstractVector{<:Real};
    bandwidth::Integer = UInt8(cld(sqrt(length(signal)), 2)),
    verbose::Bool = false,
)
    medians = Vector{Float64}(undef, Integer(cld(length(signal), bandwidth)))
    for i in eachindex(medians)
        medians[i] = median(
            signal[round(UInt, max((i - 1) * bandwidth, firstindex(signal))):round(
                UInt,
                min(i * bandwidth, lastindex(signal)),
            )],
        )
    end
    pvalue_KS = pvalue(
        ApproximateTwoSampleKSTest(
            medians[begin:round(UInt, end / 2)],
            medians[round(UInt, end / 2):end],
        ),
    )
    pvalue_JB = pvalue(JarqueBeraTest(signal; adjusted = true))
    if pvalue_KS > 0.05 || pvalue_JB > 0.05
        println("no signal detected")
    else
        z = Array{Float64}(undef, length(signal), 3)
        for i in eachindex(signal)
            z[i, 1] = if i < bandwidth
                geomean_zeros(@view(signal[begin:i]))
            else
                geomean_zeros(@view(signal[(i - bandwidth + 1):i]))
            end
            z[i, 2] = if i < bandwidth
                geovar_zeros(@view(signal[begin:i]))
            else
                geovar_zeros(@view(signal[(i - bandwidth + 1):i]))
            end
            z[i, 3] = if i < 3
                Inf
            else
                pvalue(OneSampleTTest(@view(z[2:i, 2]), mean(@view(z[2:(i - 1), 2]))))
            end
        end
        laser_start_ind = findmin(z[:, 3])[2]
        laser_start_time = time[laser_start_ind]
        q = quantile(z[laser_start_ind:end, 1], [0.05, 0.95])
        aerosol_arrival_ind = findfirst(≥(q[2]), z[:, 1])
        aerosol_arrival_time = time[aerosol_arrival_ind]
        q = quantile(z[aerosol_arrival_ind:end, 1], [0.05, 0.95])
        signal_start_ind = min(
            aerosol_arrival_ind + findfirst(<(q[2]), z[aerosol_arrival_ind:end, 1]),
            lastindex(time),
        )
        signal_end_ind = min(findlast(>(q[1]), z[:, 1]), lastindex(time))
        signal_start_time = time[signal_start_ind]
        signal_end_time = time[signal_end_ind]
        if signal_end_time - signal_start_time < 5
            if verbose === true
                println(
                    "WARNING: <5 seconds of signal initially detected after laser start. Signal times may be inaccurate or the signal has no decay.",
                )
            end
            slope = sign(
                GeochemistryTools._GLS(
                    @view(time[laser_start_ind:end]),
                    @view(signal[laser_start_ind:end]),
                    1,
                ).beta[2],
            )
            if slope > 0
                q = quantile(z[laser_start_ind:end, 1], [0.05, 0.95])
                aerosol_arrival_ind = findfirst(≥(q[1]), z[:, 1])
                aerosol_arrival_time = time[aerosol_arrival_ind]
                q = quantile(z[aerosol_arrival_ind:end, 1], [0.05, 0.95])
                signal_start_ind = min(
                    aerosol_arrival_ind + findfirst(<(q[1]), z[aerosol_arrival_ind:end, 1]),
                    lastindex(time),
                )
                signal_end_ind = min(findlast(>(q[2]), z[:, 1]), lastindex(time))
                signal_start_time = time[signal_start_ind]
                signal_end_time = time[signal_end_ind]
            end
        end
        gas_blank_end_ind = max(laser_start_ind - 2, firstindex(time))
        if verbose == true
            println(
                "gas blank: ",
                (time[begin], time[gas_blank_end_ind]),
                [1, gas_blank_end_ind],
            )
            println("laser start: ", (laser_start_time, laser_start_ind))
            println("aerosol arrival: ", (aerosol_arrival_time, aerosol_arrival_ind))
            println(
                "signal window: ",
                (signal_start_time, signal_end_time),
                (signal_start_ind, signal_end_ind),
            )
        end
        return (
            (gas_blank_end_ind, time[gas_blank_end_ind]),
            (laser_start_ind, laser_start_time),
            (aerosol_arrival_ind, aerosol_arrival_time),
            (signal_start_ind, signal_start_time),
            (signal_end_ind, signal_end_time),
        )
    end
end
