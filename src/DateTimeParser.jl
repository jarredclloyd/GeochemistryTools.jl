export automatic_datetime

"""
    automatic_datetime(datetime_string::AbstractString; [day_first::Bool=true])

    Attempt to automatically determine the `date_time_format` given only the string and
    day month order.

    - `day_first` should be either true (default) for dmy order, or false for myd order.

    Will determine if year is first by looking for a string of length 4 before the first
    delimiter. If year is first, will assume ymd order. Will not work with years only 2
    number length (e.g. 24/12/24 is ambiguous). In this case you will need to specify the
    `date_time_format` manually.

    Assumes a delimiter in the the date string, and `:` as the delimiter for time string.
"""
function automatic_datetime(datetime_string::AbstractString; day_first::Bool = true)
    if occursin(r"-", datetime_string) == true
        date_delim = '-'
    elseif occursin(r"/", datetime_string) == true
        date_delim = '/'
    end
    if occursin(r"(?i:AM|PM)", datetime_string) == true
        time_format = "H:M:S p"
    else
        time_format = "H:M:S"
    end
    if length(split(datetime_string, r"[-\/ ]")[1]) == 4
        date_format = "Y$(date_delim)m$(date_delim)d"
    elseif day_first === true
        date_format = "d$(date_delim)m$(date_delim)Y"
    elseif day_first === false
        date_format = "m$(date_delim)d$(date_delim)Y"
    end
    return DateFormat(date_format * " " * time_format)
end

function date_format_test(
    file::AbstractString;
    date_time_constructor::AbstractString,
    day_first::Bool = true,
)
    lines = split(readuntil(file, "Time "), "\n")
    date_time = rstrip(
        lines[3][(findfirst(":", lines[3])[1] + 2):(findfirst("using", lines[3])[1] - 1)],
    )
    if occursin("auto", date_time_constructor) === true
        date_time_format = automatic_datetime(date_time; day_first = day_first)
    else
        date_time_format = Dates.DateFormat(date_time_constructor)
    end
    try
        DateTime(date_time, date_time_format)
    catch err
        throw(
            ArgumentError(
                "Date format is incorrect, if specifying a custom format please see: \n \
                ? Dates.DateFormat or use automatic construction.\n If using automatic \
                construction, specify the correct format instead, \n or file a bug report \
                if it is an unambiguous string.",
            ),
        )
    end
    return date_time_format
end
