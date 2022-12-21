#= Preamble
This file contains various decay system equations for calculation of ratios and ages (U-Pb).
=#

export load_Raman, fit_base, plot_Raman_cor, plot_Raman_cor!

"""
    loadRaman(hostdir, sample; [firstrow, trailing_rows, process])

Loads Raman spectra data from txt file in hostdir.

# Description
Sample is a string representing a unique substring of the desired filename (e.g., "G17560-532nm"). 
Optional parameter `firstrow` defaults to 11 (first data row in files from [RRRUF.info](https://rruff.info/) database). 
Adjust to the required value for your input file.
Optional parameter `trailing_rows` defaults to 0, for RRRUF data change to 5. Adjust to the required value for your 
input file.
Optional parameter process defaults to true and will calculate the baseline, corrected intensities, and normalised 
intensities using the function `fit_base()`.

# Example
```julia-repl
julia> load_Raman('path/to/dir', "G17560"; firstrow = 58, trailing_rows = 0, process = true)
```
"""
function loadRaman(hostdir, sample::String;
    firstrow::Int=11,
    trailing_rows::Int=0,
    process::Bool=true)
    file = glob(sample * "*.txt", hostdir)
    data = CSV.read(file, DataFrame; header=[:wavenumber, :intensity], skipto=firstrow, types=Float64, footerskip=trailing_rows,
        ignoreemptyrows=true)
    if process == true
        fit_base(data; intensity_col=:intensity)
    elseif process == false
        return data
    end
end

"""
    fit_base(DataFrame; [intensity_col, λ])

Computes baseline and corrects intensities using the IarPLS algorithm.

# Description
Utilises pybaselines, as such requires PyCall and Conda to be installed. The GeochemistryTools.jl package will have
addressed this on package add.
`λ` (the smoothing parameter) and `intensity_col` are optional arguments which can be specified. 
The column index `(intensity_col)` parameter can be a string, integer, or symbol. It is set to :intensity by default 
which will work for data loaded using the `load_Raman()` function.
λ default is 1e5, higher numbers increase smoothing. 

Returns the original data frame with new columns (`:baseline, :corr_intensity, :norm_intensity`).

# Example
```julia-repl
julia> fit_base(df; :intensity, λ = 1e5)

```
# References
Erb, D. (2022). pybaselines: A Python library of algorithms for the baseline correction of experimental data. 
Zenodo. https://doi.org/10.5281/zenodo.7255880
"""
function fit_base(data::DataFrame; intensity_col=:intensity, λ=Nothing)
    if λ == Nothing
        output = pybaselines.whittaker.iarpls(data[!, intensity_col])
        data.baseline = output[1]
        data.corr_intensity = data[!, intensity_col] .- data[!, :baseline]
        data.norm_intensity = data[!, :corr_intensity] ./ maximum(data[!, :corr_intensity])
    elseif λ > 0
        output = pybaselines.whittaker.iarpls(data[!, intensity_col], λ)
        data.baseline = output[1]
        data.corr_intensity = data[!, intensity_col] .- data[!, :baseline]
        data.norm_intensity = data[!, :corr_intensity] ./ maximum(data[!, :corr_intensity])
    end
    return data
end

"""
    plot_raman_cor(DataFrame; [normalised, linewidth, colour])

Quickly create a lineplot of baseline corrected raman spectra. Requires use of either the load_Raman(; process = true) 
or fit_base() functions to find the correct column headers. 

# Description
Uses Makie.jl to create the plot (GLMakie is initialised on package load). Switch to CairoMakie for creating plots for 
publications. Creates a new plot and axis.
Optional parameter normalised is specified as either true or false (default is true). True will plot the normalised 
baseline-corrected intensities, whereas false will plot the baseline-corrected intensities.
Optional parameters linewidth (1) and colour ("#800080") can be set to change the linewidth and colour of the line.

# Example
```julia-repl
julia> plot_raman_cor(df; normalised = true, linewidth = 0.5, colour = :blue)
```
"""
function plot_Raman_cor(data::DataFrame;
    normalised::Bool=true,
    linewidth=1,
    colour="#800080")
    if normalised == true
        lines(data[!, :wavenumber], data[!, :norm_intensity], color=colour, linewidth=linewidth;
            axis=(; title="Normalised baseline corrected (iarPLS) Raman spectrum", xlabel="v̄ (cm⁻¹)", ylabel="Normalised intensity"))
    elseif normalised == false
        lines(data[!, :wavenumber], data[!, :corr_intensity], color=colour, linewidth=linewidth;
            axis=(; title="Baseline corrected (iarPLS) Raman spectrum", xlabel="v̄ (cm⁻¹)", ylabel="Intensity"))
    end
end

"""
    plot_raman_cor!(DataFrame; [normalised, linewidth, colour])

Create an additional lineplot of baseline corrected raman spectra on the currently active figure. 

# Description
Requires use of either the load_Raman(; process = true) or fit_base() functions to find the correct column headers. 
Uses Makie.jl to create the plot (GLMakie is initialised on package load). Switch to CairoMakie for creating plots for 
publications. Adds to the currently active plot and axis.
Optional parameter normalised is specified as either true or false (default is true). True will plot the normalised 
baseline-corrected intensities, whereas false will plot the baseline-corrected intensities.
Optional parameters linewidth (1) and colour ("#800080") can be set to change the linewidth and colour of the line.

# Example
```julia-repl
julia> plot_raman_cor!(df; normalised = true, linewidth = 0.5, colour = :blue)
```
"""
function plot_Raman_cor!(data::DataFrame;
    normalised::Bool=true,
    linewidth=1,
    colour="#800080")
    if normalised == true
        lines!(data[!, :wavenumber], data[!, :norm_intensity], color=colour, linewidth=linewidth)
    elseif normalised == false
        lines!(data[!, :wavenumber], data[!, :corr_intensity], color=colour, linewidth=linewidth)
    end
end

#= Functions in development (native Julia ports of IarPLS, ANM-IRPLS)
"""
    ANMIRPLS(ṽ, Ι)

Compute the baseline subtracted raman-shift spectra.

ṽ is a vector of the wavenumbers and Ι is a vector of intensities.
Uses the adaptive noise model based iteratively reweighted penalized least squares algorithm (ANM-IRPLS) of
Saveliev et al. 2022 (https://doi.org/10.1002/jrs.6275). 
"""
function anmirpls(ṽ, Ι)
    #1. basis construction
    fᵢ(x) = exp(-0.5 * (x - mᵢ)^2/σ^2)
    #2. regularisation

    #3. weights

    ωᵢ = 1/(10e-4 .+ Ys₁ .+ 10 .* Ys₂)^2
    #4. matrix for coefficient calculation

    #5. baseline approximation - iterative


    #6. correction of vertical displacement with noise
    m = median(Y₀ .- fit)
    bias = median((Y₀ .- fit) .* ind)
    #7. final baseline subtraction
    corrected = Y₀ .- fit - bias
end

"""
    IarPLS(v̄, I)

    Compute the baseline subtracted raman-shift spectra via IarPLS.

ṽ is a vector of the wavenumbers and Ι is a vector of intensities.
References
        ----------
        Ye, J., Tian, Z., Wei, H., & Li, Y. (2020). Baseline correction method based on improved asymmetrically 
        reweighted penalized least squares for the Raman spectrum. Applied Optics, 59(34), 10933–10943.
        https://doi.org/10.1364/AO.404863
        ---------

"""
function IarPLS(I, λ = 1e5, diff = 2, max_iter = 100, ϵ = 1e-6)
    n = length(I)
    t = 1
    Wⁱ = LinearAlgebra.I(n)
    for t ∈ 1:max_iter
        Z = (W + λ * D' * D) ^-1 * W * Y

end

function ISRU(x)
    1 / 2 * (1 - (ℯ^t * (x - 2 * σₓ₋) / σₓ₋))

end
=#