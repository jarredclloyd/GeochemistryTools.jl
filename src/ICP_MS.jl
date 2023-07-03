#=
This file contains tools for (Agilent) ICP-MS data
=#

export load_downhole, plot_dh_scatter

"""
    load_downhole(args...; kwargs...)

Load and prepare data from CSV agilent CSV files to assess downhole fractionation.

Load all relevant files in the `host_directory` using glob, concatenates the data into one table, calculates ratio between
`cps_column1` and `cps_column2`, the median-normalised ratio, and sorts the table by time.

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
    K39 -> 39, K39 -> 58) specify the desired mass by replacing ` -> ` with `_` (e.g. `K39 -> 39` to `K39_39`).

# Keywords
- `firstrow::Int`: An `Int` representing the first data row. Default value is `5`. This can be used to remove the rows
    prior to `stabletime` if set to an appropriate value (e.g, ~110 for a 30 second gas blank).
- `stabletime::Real`: The time (in seconds) where signal counts are stabilised. Default value is `32` (i.e. 32 seconds).
    This parameter is used to filter the rows when calculating the median counts per second ratio.
- `signalend::Real`: The time (in seconds) where the data should be truncated. Default value is `Inf` (will use entire
    signal beyond `stabletime`). It is recommended to set this to a value appropriate to your data (i.e. total signal
    time minus ~5 [e.g. 65]) to avoid some artefacts  that may occur near the end of an ablation.
- `trim::Bool`: Whether to trim the data to only contain values between `firstrow` and `signalend`. Default value is
    `false`. Setting to `true` it will delete all data from the first record past `signalend` to the last record.

# Example
```
julia> load_downhole("path/to/dir", "sample_x", "Rb85", "Sr87"; firstrow = 110, stabletime = 32, signalend = 65, trim = true)
MxN DataFrame

```
"""
function load_downhole(
    host_directory::AbstractString,
    sample::AbstractString,
    cps_column1::AbstractString,
    cps_column2::AbstractString;
    firstrow::Int = 5,
    stabletime::Real = 32,
    signalend::Real = Inf,
    trim::Bool = false,
)
    files = glob(sample * "*.csv", host_directory)
    data = CSV.read(
        files,
        DataFrame;
        header = 4,
        skipto = firstrow,
        types = Float64,
        footerskip = 3,
        ignoreemptyrows = true,
        normalizenames = true,
    )
    data = select!(data, r"Time", r"" * cps_column1, r"" * cps_column2)
    rename!(data, ["time", cps_column1, cps_column2])
    sort!(data, :time)
    data.ratio = data[!, cps_column1] ./ data[!, cps_column2]
    replace!(data[!, :ratio], Inf => 0)
    mdn = median(filter(x -> x .> 0, data[stabletime .< data.time .< signalend, :ratio]))
    data.ratio_mdn_norm = data[!, :ratio] ./ mdn
    if trim == true
        filter!(:time => x -> x .< signalend, data)
    end
    metadata!(data, "name", sample)
    metadata!(data, "stable time", stabletime)
    metadata!(data, "signal end", signalend)
    return data
end
