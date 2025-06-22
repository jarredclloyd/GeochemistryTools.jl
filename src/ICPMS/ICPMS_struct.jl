mutable struct LAICPMSAnalysis
    sample::Union{Nothing,AbstractString}
    material::Union{Nothing,AbstractString}
    analysis_time::Union{Nothing,DateTime}
    analysis_name::Union{Nothing,AbstractString}
    analysis_type::Union{Nothing,AbstractString}
    laser_fluence::Union{Nothing,Quantity}
    laser_repetition_rate::Union{Nothing,Quantity}
    spot_diameter::Union{Nothing,Quantity}
    laser_on::Union{Nothing,Tuple{Integer,Real}}
    stable_signal::Union{Nothing,Tuple{Integer,Real}}
    gas_blank_start::Union{Nothing,Tuple{Integer,Real}}
    gas_blank_end::Union{Nothing,Tuple{Integer,Real}}
    data::Union{Nothing,AbstractDataFrame}
    gas_blanks::Union{Nothing,Vector{Tuple{Symbol,Real,Real}}}
    signal_start::Union{Nothing,Tuple{Integer,Real}}
    signal_end::Union{Nothing,Tuple{Integer,Real}}
end

mutable struct LAICPMSSession
    laboratory::AbstractString
    session_date_start::Date
    session_date_end::Date
    ICPMS_model::AbstractString
    laser_model::AbstractString
    laser_wavelength::Quantity
    laser_pulse_width::Quantity
    analyses::Vector{LAICPMSAnalysis}
end
