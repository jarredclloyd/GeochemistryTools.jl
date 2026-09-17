# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

# GeochemistryToolsCore.jl/ReferenceValueConstants.jl
# Common constants for the GeochemistryTools ecosystem in Julia

export RADIONUCLIDES, REFERENCE_MATERIALS, GEOCHEM_REFVALS
export decay_constant, half_life, convert_ratio
export update_reference_dictionaries!, update_reference_files_remote!
export update_radionuclides!, update_geochem_refvals!, update_reference_materials!

#-----------------------------------------------------------------------------------------#
#                                 File checking
#-----------------------------------------------------------------------------------------#
# Lookup the files required for generation of the GEOCHEM_REFVALS, RADIONUCLIDES, and REFERENCE_MATERIALS constants. Will download base files from "https://github.com/jarredclloyd/GeochemistryToolsConstants/" if needed.


function _check_ref_files()
    if !ispath(_file_lookup_dir)
        mkpath(_file_lookup_dir)
    end
    if !ispath(_path_radionuclides)
        @info """"Radionuclides.yaml" missing from $_file_lookup_dir \n Downloading reference copy from github.com/jarredclloyd/GeochemistryToolsConstants"""
        Downloads.download(
            "https://raw.githubusercontent.com/jarredclloyd/GeochemistryToolsConstants/refs/heads/main/Radionuclides.yaml",
            _path_radionuclides
        )
    else
        @info "Radionuclides.yaml already exists"
    end

    if !ispath(_path_refvals)
        @info """"GeochemicalReferenceValues.yaml" missing from $_file_lookup_dir \n Downloading reference copy from github.com/jarredclloyd/GeochemistryToolsConstants"""
        Downloads.download(
            "https://raw.githubusercontent.com/jarredclloyd/GeochemistryToolsConstants/refs/heads/main/GeochemicalReferenceValues.yaml",
            _path_refvals
        )
    else
        @info "GeochemicalReferenceValues.yaml already exists"
    end

    if !ispath(_path_refmats)
        @info """"GeochemistryReferenceMaterials.yaml" missing from $_file_lookup_dir \n Downloading reference copy from github.com/jarredclloyd/GeochemistryToolsConstants"""
        Downloads.download(
            "https://raw.githubusercontent.com/jarredclloyd/GeochemistryToolsConstants/refs/heads/main/GeochemistryReferenceMaterials.yaml",
            _path_refmats
        )
    else
        @info "GeochemistryReferenceMaterials.yaml already exists"
    end
end

#-----------------------------------------------------------------------------------------#
#                                 updater functions
#-----------------------------------------------------------------------------------------#
"""
    update_radionuclides!()

    Function to update `RADIONUCLIDES` in place so a user can make alterations to $(normpath(_file_lookup_dir, "Radionuclides.yaml")) in real time.
"""
function update_radionuclides!()
    return RADIONUCLIDES[] = _construct_radionuclides(_path_radionuclides)
end

"""
    update_geochem_refvals!()

    Function to update `GEOCHEM_REFVALS` in place so a user can make alterations to $(normpath(_file_lookup_dir, "GeochemicalReferenceValues.yaml")) in real time.
"""
function update_geochem_refvals!()
    return GEOCHEM_REFVALS[] = _construct_reference_values(_path_refvals)
end

"""
    _update_reference_materials!()

    Function to update `REFERENCE_MATERIALS` in place so a user can make alterations to $(normpath(_file_lookup_dir, "GeochemistryReferenceMaterials.yaml")) in real time.
"""
function update_reference_materials!()
    return REFERENCE_MATERIALS[] = _construct_reference_materials(_path_refmats)
end

"""
    update_reference_dictionaries!()

    Function to update `RADIONUCLIDES`, `GEOCHEM_REFVALS`, and `REFERENCE_MATERIALS` in place so a user can make alterations to the files in $_file_lookup_dir in real time.
"""
function update_reference_dictionaries!()
    update_geochem_refvals!()
    update_reference_materials!()
    update_radionuclides!()
    return println("done")
end

"""
    update_reference_files_remote!()

    Function to update the local copies of "GeochemistryReferenceMaterials.yaml", "GeochemicalReferenceValues.yaml", and "Radionuclides.yaml" located in $_file_lookup_dir with the remote copies at https://github.com/jarredclloyd/GeochemistryToolsConstants

    WARNING: Will OVERWRITE existing files, will ask for confirmation.
"""
function update_reference_files_remote!()
    println("This will overwrite the following files located in $_file_lookup_dir\n
        GeochemistryReferenceMaterials.yaml\n
        GeochemicalReferenceValues.yaml\n
        Radionuclides.yaml\n
    Do you want to proceed [Y/N]")
    answer = readline()
    if occursin("y", lowercase(answer))
        _update_reference_files!()
    else
        println("aborting file update")
    end
end

function _update_reference_files!()
    if !ispath(_file_lookup_dir)
        mkpath(_file_lookup_dir)
    end
    @info "Downloading reference copy from github.com/jarredclloyd/GeochemistryToolsConstants"
    Downloads.download(
        "https://raw.githubusercontent.com/jarredclloyd/GeochemistryToolsConstants/refs/heads/main/Radionuclides.yaml",
        _path_radionuclides
    )

    @info "Downloading reference copy from github.com/jarredclloyd/GeochemistryToolsConstants"
    Downloads.download(
        "https://raw.githubusercontent.com/jarredclloyd/GeochemistryToolsConstants/refs/heads/main/GeochemicalReferenceValues.yaml",
        _path_refvals
    )

    @info "Downloading reference copy from github.com/jarredclloyd/GeochemistryToolsConstants"
    return Downloads.download(
        "https://raw.githubusercontent.com/jarredclloyd/GeochemistryToolsConstants/refs/heads/main/GeochemistryReferenceMaterials.yaml",
        _path_refmats
    )

end




#-----------------------------------------------------------------------------------------#
#                                 Half live structs and functions
#-----------------------------------------------------------------------------------------#

@enum DecayMode begin
    alpha
    beta_minus
    beta_plus
    electron_capture
    gamma
end

const DECAY_MODE_MAP = Dict(
    "alpha" => alpha,
    "α" => alpha,
    "beta_minus" => beta_minus,
    "beta_negative" => beta_minus,
    "β-" => beta_minus,
    "β⁻" => beta_minus,
    "beta_plus" => beta_plus,
    "beta_positive" => beta_plus,
    "β+" => beta_plus,
    "β⁺" => beta_plus,
    "e+" => beta_plus,
    "e⁺" => beta_plus,
    "electron_capture" => electron_capture,
    "ec" => electron_capture,
    "ϵ" => electron_capture,
    "gamma" => gamma,
    "γ" => gamma
)

Base.@kwdef struct DecayChannel <: Any
    branch::Measurement{Float64}
    mode::DecayMode
end


Base.@kwdef struct Radionuclide <: Nuclide
    radionuclide::Symbol
    references::Vector{String}
    decay_constant::Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity}
    half_life::Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity}
    decays::Base.ImmutableDict{Symbol,DecayChannel}
end

function radionuclide(;
    radionuclide::Symbol,
    references::Vector{String},
    decay_constant::Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity},
    decays::Base.ImmutableDict{Symbol,DecayChannel}
)

    half_life = Base.ImmutableDict([k => log(2) / v for (k, v) ∈ decay_constant]...)

    return Radionuclide(radionuclide, references, decay_constant, half_life, decays)
end


function Base.show(io::IO, rn::Radionuclide)
    println(io, "Radionuclide: $(rn.radionuclide)")
    println(io, "Decay constant (adopted): $(rn.decay_constant[:adopted])")
    println(io, "Half life (adopted): $(rn.half_life[:adopted])")
    println(io, "Decay mode(s): $(rn.decays)")
    return println(io, "References: \n", [ref * " \n" for ref ∈ rn.references]...)
end

"""
    _construct_radionuclides(yamlfile::String)

Internal method used to construct the `const RADIONUCLIDES` from the `Radionuclides.yaml` file that should be located at `$(normpath(_file_lookup_dir))`.

You can use this method directly to construct an alternate set of values.
It will generate a `Base.ImmutableDict::{Symbol, Radionuclide}`.
"""
function _construct_radionuclides(yamlfile::String)
    file = YAML.load_file(normpath(yamlfile); dicttype = Dict{Symbol,Any})
    return Base.ImmutableDict([
        Pair(
            key,
            radionuclide(;
                radionuclide = key,
                references = file[key][:references],
                decay_constant = Base.ImmutableDict([
                    Pair(
                        subkey,
                        measurement(
                            file[key][:decay_constant][subkey][:value],
                            file[key][:decay_constant][subkey][:uncertainty]
                        ) * uparse(
                            file[key][:decay_constant][subkey][:units];
                            unit_context = [Unitful, GeochemistryToolsCore]
                        )
                    ) for subkey ∈ keys(file[key][:decay_constant])
                ]...),
                decays = Base.ImmutableDict([
                    Pair(
                        Symbol(file[key][:decays][i][:daughter]),
                        DecayChannel(;
                            branch = measurement(
                                file[key][:decays][i][:branch][:value],
                                file[key][:decays][i][:branch][:uncertainty]
                            ),
                            mode = DECAY_MODE_MAP[file[key][:decays][i][:mode]]
                        )
                    ) for i ∈ eachindex(file[key][:decays])
                ]...)
            )
        ) for key ∈ keys(file)
    ]...)
end

"""
    RADIONUCLIDES[][:key]

    Global constant of type `Base.ImmutableDict{Symbol,GeochemistryToolsCore.Radionuclide}` that holds the `Radionuclide` structs containing some common isotopic half lives used in geochronology along with decay mode and reference information.

Access to `Radionuclide` structs within the dictionary requires suffixing `RADIONUCLIDES` with `[]` and indexing by `Symbol`.

  - e.g. RADIONUCLIDES[][:U238]

The first empty square bracket is required to enable the user to be able to adjust values in the YAML file (create new entries, use a different half life value) located at $(_file_lookup_dir)

Can be updated by explicitly calling `GeochemistryToolsCore.update_geochem_refvals!()`

Alternative half life values are available for all radioisotopes within the IsotopeTable.jl package that has been `reexported` for your use.
"""
const RADIONUCLIDES = Ref{Base.ImmutableDict{Symbol,GeochemistryToolsCore.Radionuclide}}()

"""
    decay_constant(I::Symbol; [key::Symbol=:adopted, units::String="s^-1"])



    Calculate the decay constant `λ` given a `Radionuclide` or ` Isotope` struct in `units`.

# Notes

Access `Radionuclide` via `RADIONUCLIDES[][I]` where `I` is the element symbol followed by isotope number as a `Symbol`

  - e.g.: `Symbol(U238)`, `:U238`

Provide `units` as a string compatible with Unitful.jl and `s^-1`

  - e.g.: s^-1, yr^-1, kyr^-1, Myr^-1, Gyr^-1, Maʸʳ^-1...

# Examples

```julia
julia> decay_constant(RADIONUCLIDES[][:Rb87]; units = \"Maʸʳ^-1\")
1.3972e-5 ± 4.5e-8 Ma^-1

decay_constant(isotopes[:Rb87]; units="Maʸʳ^-1")
1.3946e-5 ± 8.4e-8 Ma^-1
```
"""
function decay_constant(I::Symbol; key::Symbol = :adopted, units::String = "s^-1")
    R = RADIONUCLIDES[]
    @assert haskey(R, I) "Key ($I) not found in RADIONUCLIDES[]"
    @assert haskey(R[I][:decay_constant], key) "Key ($key) not found in decay constants for RADIONUCLIDES[][$I][:decay_constant]"
    return uconvert(
        uparse(units; unit_context = [Unitful, GeochemistryToolsCore]),
        R[I][:decay_constant][key]
    )
end

function decay_constant(h::Number; units::String = "s^-1")
    return uconvert(
        uparse(units; unit_context = [Unitful, GeochemistryToolsCore]),
        log(2) / h
    )
end

"""
    half_life(d::U where {U<:Unitful.AbstractQuantity}}; [units::String="s"])


    Calculate the `half life` given a  decay constant `λ` in `units`.

# Notes

Provide `units` as a string compatible with Unitful.jl

  - e.g.: s, yr, kyr, Myr, Gyr

# Examples

```julia
julia> half_life(measurement(1.3972e-11, 4.5e-14) * uparse(\"yr\"); units = \"Gyr\")
49.61 ± 0.16 Gyr

julia> half_life(measurement(1.3972e-11, 4.5e-14) * uparse(\"yr\"))
1.5656e18 ± 5.0e15 s
```
"""
function half_life(d::Unitful.AbstractQuantity; units::String = "s")
    return uconvert(
        uparse(units; unit_context = [Unitful, GeochemistryToolsCore]),
        log(2) / d
    )
end

#-----------------------------------------------------------------------------------------#
#                        Reference Material structs and functions
#-----------------------------------------------------------------------------------------#
# Data construction for reference materials
struct IsotopeRatio
    value::Measurement{<:Real}
    type::String
    isotope1::Symbol
    isotope2::Symbol

    function IsotopeRatio(value, type, isotope1, isotope2)
        return if isfinite(value) && ≤(value, 0)
            error("ratio cannot be ≤ 0 or non-finite, value: $value")
        else
            new(value, type, isotope1, isotope2)
        end
    end
end

"""
    isotope_ratio()

    Function used to construct IsotopeRatio `structs`,

    Accepts either positional OR keyword arguments:
    - isotope1::Union{String,Symbol},
    - isotope2::Union{String,Symbol},
    - value::Real,
    - uncertainty::Real;
    - type::Union{Nothing,String} = nothing,
"""
function isotope_ratio(
    isotope1::Union{String,Symbol},
    isotope2::Union{String,Symbol},
    value::Real,
    uncertainty::Real;
    type::Union{Nothing,String} = nothing
)
    @assert >(value, 0) "The value of an isotopic ratio cannot be ≤ 0"
    @assert isfinite(value) "A ratio that is non-finite is undefined"

    isotope1 = Symbol(isotope1)
    isotope2 = Symbol(isotope2)
    @assert in(isotope1, keys(isotopes.bysymbol)) "Isotope $isotope1 is ill-defined, key not found in `isotopes.bysymbol`"
    @assert in(isotope2, keys(isotopes.bysymbol)) "Isotope $isotope2 is ill-defined, key not found in `isotopes.bysymbol`"
    ratio = measurement(value, uncertainty)
    type = isnothing(type) ? "unspecified" : type
    return IsotopeRatio(ratio, type, isotope1, isotope2)
end

Base.@kwdef struct ReferenceMaterial
    ID::Symbol
    IGSN::Union{Nothing,String}
    name::String
    description::String
    material::String
    aliases::Union{Nothing,String,Vector{String}} = nothing
    notes::Union{Nothing,String,Vector{String}} = nothing
    references::Union{Nothing,String,Vector{String}} = nothing
    elemental_concentrations::Union{
        Nothing,
        Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity}
    } = nothing
    isotopic_abundances::Union{Nothing,Base.ImmutableDict{Symbol,Measurement{Float64}}} =
        nothing
    isotopic_ratios::Union{Nothing,Base.ImmutableDict{Symbol,IsotopeRatio}} = nothing
    age_data::Union{Nothing,Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity}} = nothing
end

function Base.show(io::IO, RM::ReferenceMaterial)
    println("Key: $(RM.ID)")
    if !isnothing(RM.IGSN)
        println("IGSN: https:://doi.org/$(RM.IGSN)")
    end
    println("Name: $(RM.name)")
    println("Material: $(RM.material)")
    println(RM.description)
    if !isnothing(RM.aliases)
        println("Aliases: $(join(RM.aliases, ", "))")
    end
    if !isnothing(RM.notes)
        println("Notes: $(RM.notes)")
    end
    if !isnothing(RM.references)
        println("References: $(join(RM.references, "; "))")
    end
    if !isnothing(RM.elemental_concentrations)
        println("Available elemental concentrations: $(join(keys(RM.elemental_concentrations), ", "))",)
    end
    if !isnothing(RM.isotopic_abundances)
        println("Specified isotopic abundances: $(join(keys(RM.isotopic_abundances), ", "))",)
    end
    if !isnothing(RM.isotopic_ratios)
        println("Available isotopic ratios: $(join(keys(RM.isotopic_ratios), ", "))")
    end
    if !isnothing(RM.age_data)
        println("Available age determinations: $(join(keys(RM.age_data), ", "))")
    end
end

"""
    reference_material()

    Functions used to construct ReferenceMaterial `structs`.

    Accepts either positional OR keyword arguments:
    - ID::Symbol
    - name::String
    - material::String
    - aliases::Union{Nothing,String,Vector{String}} = nothing
    - references::Union{Nothing,String,Vector{String}} = nothing
    - elemental_concentrations::Union{Nothing,Base.ImmutableDict{Symbol,Tuple{Float64,Float64}},} = nothing
    - isotopic_abundances::Union{Nothing,Base.ImmutableDict{Symbol,Tuple{Float64,Float64}}} = nothing
    - isotopic_ratios::Union{Nothing,Base.ImmutableDict{Symbol,IsotopeRatio}} = nothing
    - age_data::Union{Nothing,Base.ImmutableDict{Symbol,Tuple{Float64,Float64}}} = nothing
"""
function reference_material(
    ID::Symbol,
    IGSN::Union{Nothing,String},
    name::String,
    description::String,
    material::String,
    aliases::Union{Nothing,String,Vector{String}} = nothing,
    notes::Union{Nothing,String,Vector{String}} = nothing,
    references::Union{Nothing,String,Vector{String}} = nothing,
    elemental_concentrations::Union{
        Nothing,
        Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity}
    } = nothing,
    isotopic_abundances::Union{Nothing,Base.ImmutableDict{Symbol,<:Measurement}} = nothing,
    isotopic_ratios::Union{Nothing,Base.ImmutableDict{Symbol,IsotopeRatio}} = nothing,
    age_data::Union{Nothing,Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity}} = nothing
)
    return ReferenceMaterial(;
        ID = ID,
        IGSN = IGSN,
        name = name,
        description = description,
        material = material,
        aliases = aliases,
        notes = notes,
        references = references,
        elemental_concentrations = elemental_concentrations,
        isotopic_abundances = isotopic_abundances,
        isotopic_ratios = isotopic_ratios,
        age_data = age_data
    )
end

"""
    _construct_reference_values(yamlfile::String)

Internal method used to construct the `const REFERENCE_MATERIALS` from the `GeochemistryReferenceMaterials.yaml` file that should be located at `$(normpath(_file_lookup_dir))`.

You can use this method directly to construct an alternate set of values.
It will generate a `Base.ImmutableDict::{Symbol, ReferenceMaterial}`.
"""
function _construct_reference_materials(yamlfile::String)
    file = YAML.load_file(normpath(yamlfile); dicttype = Dict{Symbol,Any})
    return Base.ImmutableDict([
        Pair(
            Symbol(file[key][:refmatID]),
            ReferenceMaterial(;
                ID = Symbol(file[key][:refmatID]),
                IGSN = isnothing(file[key][:IGSN]) ? file[key][:IGSN] : nothing,
                name = file[key][:name],
                description = file[key][:description],
                material = file[key][:material],
                aliases = file[key][:aliases],
                notes = file[key][:notes],
                references = file[key][:references],
                elemental_concentrations = if isnothing(file[key][:elemental_concentrations],)
                    nothing
                else
                    Base.ImmutableDict([
                        Pair(
                            subkey,
                            measurement(
                                file[key][:elemental_concentrations][subkey][:value],
                                file[key][:elemental_concentrations][subkey][:uncertainty]
                            ) * uparse(
                                file[key][:elemental_concentrations][subkey][:units];
                                unit_context = [Unitful, GeochemistryToolsCore]
                            )
                        ) for subkey ∈ keys(file[key][:elemental_concentrations])
                    ]...,)
                end,
                isotopic_abundances = if isnothing(file[key][:isotopic_abundances])
                    nothing
                else
                    Base.ImmutableDict([
                        Pair(
                            subkey,
                            measurement(
                                file[key][:isotopic_abundances][subkey][:value],
                                file[key][:isotopic_abundances][subkey][:uncertainty]
                            )
                        ) for subkey ∈ keys(file[key][:isotopic_abundances])
                    ]...,)
                end,
                isotopic_ratios = if isnothing(file[key][:isotopic_ratios])
                    nothing
                else
                    Base.ImmutableDict([
                        Pair(
                            subkey,
                            isotope_ratio(
                                file[key][:isotopic_ratios][subkey][:isotope1],
                                file[key][:isotopic_ratios][subkey][:isotope2],
                                file[key][:isotopic_ratios][subkey][:value],
                                file[key][:isotopic_ratios][subkey][:uncertainty];
                                type = file[key][:isotopic_ratios][subkey][:ratio_type]
                            )
                        ) for subkey ∈ keys(file[key][:isotopic_ratios])
                    ]...,)
                end,
                age_data = if isnothing(file[key][:age_data])
                    nothing
                else
                    Base.ImmutableDict([
                        Pair(
                            subkey,
                            measurement(
                                file[key][:age_data][subkey][:value],
                                file[key][:age_data][subkey][:uncertainty]
                            ) * uparse(
                                file[key][:age_data][subkey][:units];
                                unit_context = [Unitful, GeochemistryToolsCore]
                            )
                        ) for subkey ∈ keys(file[key][:age_data])
                    ]...,)
                end
            )
        ) for key ∈ keys(file)
    ]...,)
end

"""
    REFERENCE_MATERIALS[][:key]

    Global constant of type `Base.ImmutableDict{Symbol,GeochemistryToolsCore.ReferenceMaterial}` that holds the `ReferenceMaterial` structs containing information and values for some reference materials used in geochemistry.

Access to `ReferenceMaterial` structs within the dictionary requires suffixing `REFERENCE_MATERIALS` with `[]` and indexing by `Symbol`.

  - e.g. REFERENCE_MATERIALS[][:NIST610]

The first empty square bracket is required to enable the user to be able to adjust values in the YAML file (create new entries, use a different set of values for a reference material) located at $(_file_lookup_dir)

Can be updated by explicitly calling `GeochemistryToolsCore.update_geochem_refvals!()`
"""
const REFERENCE_MATERIALS =
    Ref{Base.ImmutableDict{Symbol,GeochemistryToolsCore.ReferenceMaterial}}()

"""
    convert_ratio(IR::IsotopeRatio, [frac_isotope1=nothing, frac_isotope2=nothing])

Convert from and atomic (molar) isotope ratio to a mass isotope ratio and vice versa.

This version of the function makes the conversion by using an intermediary to calculate as a fraction of total element.

If no value (of type `Float64`, `Rational`, or `Measurement(Float64)`) is specified for `frac_isotope1` or `frac_isotope2` the natural isotopic abundance will be used in place.
"""
function convert_ratio(
    IR::IsotopeRatio,
    frac_isotope1::Union{Nothing,Real,Measurement{Float64}} = nothing,
    frac_isotope2::Union{Nothing,Real,Measurement{Float64}} = nothing
)
    if occursin("unspecified", IR.type)
        @warn "ratio type unspecified, no conversion performed"
    else
        if isnothing(frac_isotope1)
            frac_isotope1 = isotopes[IR.isotope1].abundance / 100
            @info "assuming natural abundance for $(IR.isotope1) as no fraction specified"
        end
        if isnothing(frac_isotope2)
            frac_isotope2 = isotopes[IR.isotope2].abundance / 100
            @info "assuming natural abundance for $(IR.isotope2) as no fraction specified"
        end
        @assert (0 < frac_isotope1 ≤ 1) "frac_isotope1 must meet the condition 0 < x ≤ 1, got value: $frac_isotope1"
        @assert (0 < frac_isotope2 ≤ 1) "frac_isotope2 must meet the condition 0 < x ≤ 1, got value: $frac_isotope2"

        element1_mass = measurement(
            Mendeleev.elements[isotopes[IR.isotope1].atomic_number].atomic_weight.val,
            Mendeleev.elements[isotopes[IR.isotope1].atomic_number].atomic_weight_uncertainty
        )
        element2_mass = measurement(
            Mendeleev.elements[isotopes[IR.isotope2].atomic_number].atomic_weight.val,
            Mendeleev.elements[isotopes[IR.isotope2].atomic_number].atomic_weight_uncertainty
        )
        if occursin("mass", IR.type)
            return IsotopeRatio(
                (IR.value / frac_isotope1 / element1_mass * frac_isotope1) /
                (1 / frac_isotope2 / element2_mass * frac_isotope2),
                "atomic",
                IR.isotope1,
                IR.isotope2
            )
        elseif occursin("atomic", IR.type)
            IsotopeRatio(
                (IR.value / frac_isotope1 * element1_mass * frac_isotope1) /
                (1 / frac_isotope2 * element2_mass * frac_isotope2),
                "mass",
                IR.isotope1,
                IR.isotope2
            )
        end
    end
end

"""
    convert_ratio(IR::Union{Real, Measurement{Float64}},
    type::String,
    isotope1::Symbol,
    isotope2::Symbol,
    frac_isotope1::Union{Nothing,Float64,Rational,Measurement{Float64}} = nothing,
    frac_isotope2::Union{Nothing,Float64,Rational,Measurement{Float64}} = nothing,)

Convert from and atomic (molar) isotope ratio to a mass isotope ratio and vice versa.

This version of the function makes the conversion by using an intermediary to calculate as a fraction of total element.

If no value (of type `Float64`, `Rational`, or `Measurement(Float64)`) is specified for `frac_isotope1` or `frac_isotope2` the natural isotopic abundance will be used in place.
"""
function convert_ratio(
    IR::Union{Real,Measurement{Float64}},
    type::String,
    isotope1::Symbol,
    isotope2::Symbol,
    frac_isotope1::Union{Nothing,Real,Measurement{Float64}} = nothing,
    frac_isotope2::Union{Nothing,Real,Measurement{Float64}} = nothing
)
    @assert (occursin("atomic", type) || occursin("mass", type)) "ratio type should be either `atomic`` or `mass`"
    @assert (in(isotope1, keys(isotopes.bysymbol))) "$isotope1 not in keys(isotopes.bysymbol)"
    @assert (in(isotope2, keys(isotopes.bysymbol))) "$isotope2 not in keys(isotopes.bysymbol)"
    if isnothing(frac_isotope1)
        frac_isotope1 = isotopes[isotope1].abundance / 100
        @info "assuming natural abundance for $(isotope1) as no fraction specified"
    end
    if isnothing(frac_isotope2)
        frac_isotope2 = isotopes[isotope2].abundance / 100
        @info "assuming natural abundance for $(isotope2) as no fraction specified"
    end
    @assert (0 < frac_isotope1 ≤ 1) "frac_isotope1 must meet the condition 0 < x ≤ 1, got value: $frac_isotope1"
    @assert (0 < frac_isotope2 ≤ 1) "frac_isotope2 must meet the condition 0 < x ≤ 1, got value: $frac_isotope2"

    element1_mass = measurement(
        Mendeleev.elements[isotopes[isotope1].atomic_number].atomic_weight.val,
        Mendeleev.elements[isotopes[isotope1].atomic_number].atomic_weight_uncertainty
    )
    element2_mass = measurement(
        Mendeleev.elements[isotopes[isotope2].atomic_number].atomic_weight.val,
        Mendeleev.elements[isotopes[isotope2].atomic_number].atomic_weight_uncertainty
    )
    if occursin("mass", type)
        return IsotopeRatio(
            (IR / frac_isotope1 / element1_mass * frac_isotope1) /
            (1 / frac_isotope2 / element2_mass * frac_isotope2),
            "atomic",
            isotope1,
            isotope2
        )
    elseif occursin("atomic", type)
        IsotopeRatio(
            (IR / frac_isotope1 * element1_mass * frac_isotope1) /
            (1 / frac_isotope2 * element2_mass * frac_isotope2),
            "mass",
            isotope1,
            isotope2
        )
    end
end

#-----------------------------------------------------------------------------------------#
#                           Reference value structs and functions
#-----------------------------------------------------------------------------------------#
Base.@kwdef struct ReferenceValue
    ID::Symbol
    name::String
    class::String
    reference::String
    elemental_concentrations::Union{
        Nothing,
        Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity}
    } = nothing
end

function Base.show(io::IO, RV::ReferenceValue)
    println("Key: $(RV.ID)")
    println("Name: $(RV.name)")
    println("Class: $(RV.class)")
    println("Reference: $(RV.reference)")
    if !isnothing(RV.elemental_concentrations)
        println("Available elemental concentrations: $(join(keys(RV.elemental_concentrations), ", "))",)
    end
end

"""
    reference_value()

    Functions used to construct ReferenceValue `structs`.

    Accepts either positional OR keyword arguments:
    - ID::Symbol
    - name::String
    - class::String
    - references::String
    - elemental_concentrations::Union{Nothing,Base.ImmutableDict{Symbol,Tuple{Float64,Float64}},} = nothing
"""
function reference_value(
    ID::Symbol,
    name::String,
    class::String,
    reference::String,
    elemental_concentrations::Union{
        Nothing,
        Base.ImmutableDict{Symbol,<:Unitful.AbstractQuantity}
    } = nothing
)
    return ReferenceValue(;
        ID = ID,
        name = name,
        class = class,
        reference = reference,
        elemental_concentrations = elemental_concentrations
    )
end

"""
    _construct_reference_values(yamlfile::String)

Internal method used to construct the `const GEOCHEM_REFVALS` from the `GeochemicalReferenceValues.yaml` file that should be located at `$(normpath(_file_lookup_dir))`.

You can use this method directly to construct an alternate set of values.
It will generate a `Base.ImmutableDict::{Symbol, ReferenceValue}`.
"""
function _construct_reference_values(yamlfile::String)
    file = YAML.load_file(normpath(yamlfile); dicttype = Dict{Symbol,Any})
    return Base.ImmutableDict([
        Pair(
            Symbol(file[key][:refvalID]),
            ReferenceValue(;
                ID = Symbol(file[key][:refvalID]),
                name = file[key][:name],
                class = file[key][:class],
                reference = file[key][:reference],
                elemental_concentrations = if isnothing(file[key][:elemental_concentrations],)
                    nothing
                else
                    Base.ImmutableDict([
                        Pair(
                            subkey,
                            measurement(
                                file[key][:elemental_concentrations][subkey][:value],
                                file[key][:elemental_concentrations][subkey][:uncertainty]
                            ) * uparse(
                                file[key][:elemental_concentrations][subkey][:units];
                                unit_context = [Unitful, GeochemistryToolsCore]
                            )
                        ) for subkey ∈ keys(file[key][:elemental_concentrations])
                    ]...,)
                end
            )
        ) for key ∈ keys(file)
    ]...,)
end

"""
    GEOCHEM_REFVALS[][:key]

    Global constant of type `Base.ImmutableDict{Symbol,GeochemistryToolsCore.ReferenceValue}` that holds the `ReferenceValue` structs containing some common isotopic half lives used in geochronology along with decay mode and reference information.

Access to `ReferenceValue` structs within the dictionary requires suffixing `GEOCHEM_REFVALS` with `[]` and indexing by `Symbol`.

  - e.g. GEOCHEM_REFVALS[][:CI_CHONDRITE_PO2014

The first empty square bracket is required to enable the user to be able to adjust values in the YAML file (create new entries, use a different half life value) located at $(_file_lookup_dir)

Can be updated by explicitly calling `GeochemistryToolsCore.update_geochem_refvals!()`
"""
const GEOCHEM_REFVALS =
    Ref{Base.ImmutableDict{Symbol,GeochemistryToolsCore.ReferenceValue}}()
