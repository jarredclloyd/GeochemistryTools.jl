#=
Functions for barycentric coordinate conversions
=#

export cartesian, ternary

"""
    cartesian(a, b, c)

Converts ternary coordinates (a, b, c) to cartesian coordinates (x, y).

```
julia> cartesian(-1.1547005383792517, -0.15470053837925168, 2.3094010767585034)
(1, 2)
```
"""
function cartesian(a::Real, b::Real, c::Real)
    x = 0.5 * (2b + c) / (a + b + c)
    y = √3 / 2 * (c / (a + b + c))
    return x, y
end

"""
    ternary(x, y)

Converts cartesian coordinates (x, y) to ternary coordinates (a, b, c).

```
julia> ternary(1,2)
(-1.1547005383792517, -0.15470053837925168, 2.3094010767585034)
```
"""
function ternary(x::Real, y::Real)
    c = (2 * y) / √3
    b = x - c / 2
    a = 1 - b - c
    return a, b, c
end

"""
    ternary(a, b, c, d)

Converts tetrahedral coordinates (a, b, c, d) to ternary coordinates (a, b, c).

```
julia> ternary(1, 0, 0, 1)
(0.5, 1.1547005383792515, 0.8164965809277259)
```
"""
function ternary(a::Real, b::Real, c::Real, d::Real)
    x = (c + 1 - b)/2
    y = √3 / 2 * a + √3 / 6 * d
    z = √6 / 3 * d
    return (x, y, z)
end
