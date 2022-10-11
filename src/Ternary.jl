#=
Functions for plotting ternary plots
=#

export ternary2cartesian, cartesian2ternary

"""
    ternary2cartesian(x, y)

Computes the cartesian coordinates (x, y) for a given set of ternary coordinates (a, b, c).

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

Computes the ternary coordinates (a, b, c) for a given set of cartesian coordinates (x, y).

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



