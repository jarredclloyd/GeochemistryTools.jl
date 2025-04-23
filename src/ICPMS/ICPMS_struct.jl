struct LAICPMSAnalysis
    session_date_start::Date
    material::AbstractString
    sample::AbstractString
    analysis::AbstractString
    analysis_time::Time
    laser_fluence::AbstractFloat
    laser_repetition_rate::Real
    laser_on::AbstractFloat
    laser_off::AbstractFloat
    spot_diameter::Integer
    gas_blank_start::Real
    gas_blank_end::Real
    gas_blanks::Vector{Real}
    signal_start::Tuple{Real, Integer}
    signal_end::Tuple{Real, Integer}
    data::DataFrame
end

struct LAICPMSSession
    laboratory::AbstractString
    session_date_start::Date
    session_date_end::Date
    ICPMS_model::AbstractString
    laser_model::AbstractString
    laser_wavelength::Real
    laser_pulse_width::Real
end
