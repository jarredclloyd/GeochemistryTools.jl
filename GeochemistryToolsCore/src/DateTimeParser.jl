# SPDX-FileCopyrightText: Copyright © 2026 Jarred C. Lloyd
# SPDX-FileContributor: Jarred C. Lloyd <4o8pg2v9h@mozmail.com>
#
# SPDX-License-Identifier: MPL-2.0

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

function _date_format_test(date_time::AbstractString, date_time_format::DateFormat)
    test_passed::Bool=false
    try
        DateTime(date_time, date_time_format)
        test_passed = true
    catch err
    end
    return test_passed
end
