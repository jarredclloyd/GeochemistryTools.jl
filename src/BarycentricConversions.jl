#=
Functions for barycentric coordinate conversions
=#

export ternary2cartesian, cartesian2ternary, tetradhedral2ternary

"""
    ternary2cartesian(x, y)

Converts ternary coordinates (a, b, c) to cartesian coordinates (x, y).

```
julia> cartesian2ternary(-1.1547005383792517, -0.15470053837925168, 2.3094010767585034)
(1, 2)
```
"""
function ternary2cartesian(a, b, c)
    x = 0.5 * (2b + c) / (a + b + c)
    y = √3 / 2 * (c / (a + b + c))
    return x, y
end

"""
    cartesian2ternary(x, y)

Converts cartesian coordinates (x, y) to ternary coordinates (a, b, c).

```
julia> cartesian2ternary(1,2)
(-1.1547005383792517, -0.15470053837925168, 2.3094010767585034)
```
"""
function cartesian2ternary(x, y)
    c = (2 * y) / √3
    b = x - c / 2
    a = 1 - b - c
    return a, b, c
end


function tetradhedral2ternary(a, b, c, d)
    x = (c + 1 - b)/2
    y = √3 / 2 * a + √3 / 6 * d
    z = √6 / 3 * d
    return (x, y, z)
end
