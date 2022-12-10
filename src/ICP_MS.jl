#=
This file contains tools for (Agilent) ICP-MS data
=#

export load_downhole, plot_dh_scatter

"""
    load_downhole(hostdir, sample::String, CPS_col1::String, CPS_col2::String, [firstrow::Int])

Load and prepare data from CSV agilent CSV files to assess downhole fractionation.
Works for vectors and other arrays if column is specified.

# Examples
```
julia> load_downhole("path/to/dir", "sample_x", "Rb85", "Sr87", firstrow = 110)
MxN DataFrame

```
"""
function load_downhole(hostdir, sample::String, CPS_col1::String, CPS_col2::String, firstrow=5::Int)
    files = glob(sample * "*.csv", hostdir)
    data = CSV.read(files, DataFrame; header=4, skipto=firstrow, types=Float64, footerskip=3, ignoreemptyrows=true, normalizenames=true)
    data = select!(data, r"Time", r"" * CPS_col1, r"" * CPS_col2)
    rename!(data, ["Time", CPS_col1, CPS_col2])
    sort!(data, :Time)
    transform!(data, [CPS_col1, CPS_col2] => ByRow((x, y) -> /(x, y)) => :ratio)
    replace!(data[!, :ratio], Inf => 0)
    mdn = median(filter(x -> x .> 0, data[!, :ratio]))
    transform!(data, :ratio => (x -> x / mdn) => :ratio_mdn_norm)
end


"""
    plot_dh_scatter(data, normalised=true, pointsize=5, colour="#800080")

Create a basic scatter plot of the downhole fractionation pattern for a dataset loaded via "load_downhole".
Can specify "normalised", size of maker, and colour. Only requires the dataframe name, and will default to the normalised plot.

# Examples
```
julia>  plot_dh_scatter(df)

```
"""
function plot_dh_scatter(data::DataFrame, normalised=true::Bool, pointsize=5, colour="#800080")
    if normalised == true
        scatter(data[!, :Time], data[!, :ratio_mdn_norm], color=colour, markersize=pointsize;
            axis=(; title="Median normalised downhole fractionation", xlabel="Time (s)", ylabel="(CPS/CPS)/median"))
    elseif normalised == false
        scatter(data[!, :Time], data[!, :ratio], color=colour, markersize=pointsize;
            axis=(; title="Downhole fractionation", xlabel="Time (s)", ylabel="CPS/CPS"))
    end
end