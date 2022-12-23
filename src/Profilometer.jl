#= 
Functions for laser profilometry (Olympus OLS5000)
=#

export loadProfilometer

"""
    loadProfilometer(hostdir::string, sample::string; pixelwidth::Int=1024, headerrow::Int=19, datarow::Int=20)

    Loads Olympus OLS5000 profilometer CSV data (exported as area)
"""
function loadProfilometer(hostdir, sample; headerrow::Int=19, datarow::Int=20)
    file = glob(sample * "*.csv", hostdir)
    df = CSV.read(file, DataFrame; header=false, skipto=5, limit=2, select=[1, 2])
    yresolution = df[1, 2]
    zresolution = df[2, 2]
    df = CSV.read(file, DataFrame; header=headerrow, skipto=datarow, select=(index, name) -> occursin("POS", 
        String(name)) == true || occursin("DataLine", String(name)) == true)
    pixelwidth = ncol(df)
    df = stack(df, 2:pixelwidth)
    df.POS = replace.(df[!, 2], "POS = " => "")
    df.y = 0 .+ parse.(Int, df.POS) .* yresolution
    select!(df, :DataLine => :x, :y, :value => :z)
    df.z = (df[!, :z]) .- mode(df[!, :z])
    metadata!(df, "xyresolution", yresolution); metadata!(df, "zresolution", zresolution);
    return df
end