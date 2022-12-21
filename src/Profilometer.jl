#= 
Functions for laser profilometry (Olympus OLS5000)
=#
"""
    loadProfilometer(hostdir::String, sample::String; pixelwidth::Int=1024, headerrow::Int=19, datarow::Int=20)

    Loads Olympus OLS5000 profilometer CSV data (exported as area)
"""
function loadProfilometer(hostdir, sample; pixelwidth = 1024, headerrow = 19, datarow = 20)
    file = glob(sample * "*.csv", hostdir)
    df = CSV.read(file, DataFrame; header = headerrow, skipto = datarow, drop = (index, name) -> index > pixelwidth)
    df = stack(df, 2:pixelwidth)
    df.POS = replace.(df[!, 2], "POS = " => "")
    df.y = 0 .+ parse.(Int, df.POS) .* 0.25
    select!(df, :DataLine => :x, :y, :value => :z)
    df.z = (df[!, :z]) .- maximum(df[!, :z])
    
    return df
end