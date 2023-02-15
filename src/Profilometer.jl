#= 
Functions for laser profilometry (Olympus OLS5000)
=#

export loadProfilometer

"""
    loadProfilometer(hostdir::string, sample::string; headerrow::Int=19, datarow::Int=20)

Loads Olympus OLS5000 profilometer CSV data (exported as height and area).

# Description
Loads CSV data outputs (CSV image) from Olympus profilometer software. Translates "POS = X" to appropriate coordinate
based on the resolution fields in the file. It adjusts the height to be zeroed around the median value, which should
be equal to the normal surface height prior to ablation.

# Example
```julia-repl
julia> loadProfilometer(path/to/host/dir/, "Sample_X"; headerrow = 19, datarow = 20)
MxN DataFrame
```

"""
function loadProfilometer(hostdir, sample; headerrow::Int=19, datarow::Int=20)
    file = glob(sample * "*.csv", hostdir)
    df = CSV.read(file, DataFrame; header=false, skipto=4, limit=3, select=[1, 2])
    xresolution = df[1, 2]
    yresolution = df[2, 2]
    zresolution = df[3, 2]
    df = CSV.read(file, DataFrame; header=headerrow, skipto=datarow, select=(index, name) -> occursin("POS",
        String(name)) == true || occursin("DataLine", String(name)) == true)
    pixelwidth = ncol(df)
    df = stack(df, 2:pixelwidth)
    df.POS = replace.(df[!, 2], "POS = " => "")
    df.y = 0 .+ parse.(Int, df.POS) .* xresolution
    select!(df, :DataLine => :x, :y, :value => :z)
    df.z = (df[!, :z]) .- mode(df[!, :z])
    metadata!(df, "xresolution", xresolution)
    metadata!(df, "yresolution", yresolution)
    metadata!(df, "zresolution", zresolution)
    return df
end