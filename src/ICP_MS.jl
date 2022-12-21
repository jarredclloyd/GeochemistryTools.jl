#=
This file contains tools for (Agilent) ICP-MS data
=#

export loadDownhole, plot_dh_scatter

"""
    loadDownhole(hostdir, sample, CPS_col1, CPS_col2, [firstrow, stable_time])

Load and prepare data from CSV agilent CSV files to assess downhole fractionation.

# Description
Load all relevant files in the hostdir using glob (e.g. NIST610 will just load `NIST610*.csv` files, whereas NIST 
will load any file where the name contains `NIST*.csv`). 

`CPS_cols` should be specified as element mass (e.g., Rb85). It matches using regular expression and normalises the 
header names. If duplicate isotope values are present (due to mass shift checks, e.g. K39 -> 39, K39 -> 58) specify the 
desired mass by replacing ` -> ` with `_` (e.g. `K39 -> 39` to `K39_39`).

Concatenates the data into one table, calculates ratio between `CPS_col1` and `CPS_col2`, the median-normalised ratio, 
and sorts the table by time. 

It is recommended to specify the `stable_time` parameter to the point where the signal counts are stabilised. It is set by 
default to be 32 (i.e. 32 seconds). This parameter is used to filter the rows to calculate the median CPS ratio. If you 
want to remove the rows prior to the stable signal time, adjust `firstrow` to equal the row equivalent to the stable_time
(usually about 110 for a 30 second gas blank).

# Example
```
julia> loadDownhole("path/to/dir", "sample_x", "Rb85", "Sr87"; firstrow = 5, stable_row = 110)
MxN DataFrame

```
"""
function loadDownhole(hostdir, sample::String, CPS_col1::String, CPS_col2::String;
    firstrow::Int = 5,
    stable_time::Number = 32)
    files = glob(sample * "*.csv", hostdir)
    data = CSV.read(files, DataFrame; header = 4, skipto = firstrow, types = Float64, footerskip = 3, 
        ignoreemptyrows = true, normalizenames = true)
    data = select!(data, r"Time", r"" * CPS_col1, r"" * CPS_col2)
    rename!(data, ["time", CPS_col1, CPS_col2])
    sort!(data, :time)
    data.ratio =  data[!, CPS_col1] ./ data[!, CPS_col2]
    replace!(data[!, :ratio], Inf => 0)
    mdn = median(filter(x -> x .> 0, data[data.time .> stable_time, :ratio]))
    data.ratio_mdn_norm = data[!,:ratio] ./ mdn
    metadata!(data, "name", sample); metadata!(data,"stable_time", stable_time);
    return data
end

"""
    plot_dh_scatter(DataFrame; [stable_time, laser_start, normalised, marker_size, marker_colour, line_colour])

# Description
Create a basic scatter plot of the downhole fractionation pattern for a dataset loaded via `load_downhole`.

Can specify `name`, `stable_time`, `laser_start`, `normalised`, `marker_size`, `marker_colour`, and `line_colour`. 
Only requires the DataFrame name, and will default to the median-normalised plot. 
Name and `stable_time` will be taken from the metadata of the input DataFrame if loaded using `load_downhole()`. 
`laser_start` will default to `stable_time - 2` but can be adjusted to whatever value you want the plot to start at.

# Example
```
julia>  plot_dh_scatter(df)

```
"""
function plot_dh_scatter(data::DataFrame; 
    name = metadata(data, "name"),
    stable_time = metadata(data, "stable_time"),
    laser_start = stable_time - 2,
    normalised::Bool = true,
    marker_size = 5,
    marker_colour = "#800080", 
    line_colour = :orange)
    fig1 = Figure(resolution=(1920,1080))
    if normalised == true
        fit = Polynomials.fit(data[data.time .> stable_time, :time], data[data.time .> stable_time, :ratio_mdn_norm], 3)
        ax1 = Axis(fig1[1, 1]; title = "Median normalised downhole fractionation — " * name, xlabel = "Time (s)",
            ylabel = "(CPS/CPS)/median")
        scatter!(ax1, data[data.time .> laser_start, :time], data[data.time .> laser_start, :ratio_mdn_norm],
         color = marker_colour, markersize = marker_size)
        lines!(ax1, data[data.time .> stable_time, :time], coeffs(fit)[1] .+ coeffs(fit)[2] .* 
            data[data.time .> stable_time, :time] .+ coeffs(fit)[3] .* data[data.time .> stable_time, :time] .^ 2 .+
            coeffs(fit)[4] .* data[data.time .> stable_time, :time] .^ 3, color = line_colour)
    elseif normalised == false
        fit = Polynomials.fit(data[data.time .> stable_time, :time], data[data.time .> stable_time, :ratio], 3)
        ax1 = Axis(fig1[1, 1]; title = "Downhole fractionation — " * name, xlabel ="Time (s)", ylabel="CPS/CPS")
        scatter(ax1, data[data.time .> laser_start, :time], data[data.time .> laser_start, :ratio],
         color = marker_colour, markersize = marker_size)
        lines!(ax1, data[data.time .> stable_time, :time], coeffs(fit)[1] .+ coeffs(fit)[2] .* 
            data[data.time .> stable_time, :time] .+ coeffs(fit)[3] .* data[data.time .> stable_time, :time] .^ 2 .+
            coeffs(fit)[4] .* data[data.time .> stable_time, :time] .^ 3, color = line_colour)
    end
    return fig1
end