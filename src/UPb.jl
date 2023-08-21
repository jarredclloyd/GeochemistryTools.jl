#= Preamble
This file contains various decay system equations for calculation of ratios and ages (U-Pb).
=#

export age2ratioPb207U235, age2ratioPb206U238, age2ratioPb207Pb206, ratio2agePb207U235, ratio2agePb206U238,
    ratio2agePb207Pb206, calcaitchisonTW

#Function calls
"""
    age2ratioPb206U238(age)

Compute the ²⁰⁶Pb/²³⁸U ratio for a given age.

Works for vectors and other arrays if column is specified.

# Examples
```julia-repl
julia> age2ratioPb206U238(1000)
 0.16780392747297124

julia> a = [1000; 2000; 3000];
julia> age2ratioPb206U238(a)
 0.16780392747297124
 0.3637660130212965
 0.592611306160425
```
"""
function age2ratioPb206U238(age)
    if isa(age, Array)
        ratio = zeros(size(age))
        nages = length(age)
        @simd for i in 1:nages
            ratio[i] = ratioPb206U238(age[i])
        end
        return ratio
    else
        ratio = ratioPb206U238(age)
    end
end

"""
    age2ratioPb207U235(age)

Compute the ²⁰⁷Pb/²³⁵U ratio for a given age.

Works for vectors and other arrays if column is specified.

# Examples
```julia-repl
julia> age2ratioPb207U235(1000)
 1.6774102427622641

julia> a = [1000; 2000; 3000]
julia> age2ratioPb207U235(a)
  1.6774102427622641
  6.168525608048285
 18.193083888492072
```
"""
function age2ratioPb207U235(age)
    if isa(age, Array)
        ratio = zeros(size(age))
        nages = length(age)
        @simd for i in 1:nages
            ratio[i] = ratioPb207U235(age[i])
        end
        return ratio
    else
        ratio = ratioPb207U235(age)
    end
end

"""
    age2ratioPb207Pb206(age)

Compute the ²⁰⁷Pb/²⁰⁶Pb ratio for a given age.

Works for vectors and other arrays if column is specified.

# Examples
```julia-repl
julia> age2ratioPb207Pb206(1000)
  0.07253226274578425

julia> a = [1000; 2000; 3000]
julia> age2ratioPb207Pb206(a)
 0.07253226274578425
 0.1230419792914581
 0.22275653097919715
```
"""
function age2ratioPb207Pb206(age)
    if isa(age, Array)
        ratio = zeros(size(age))
        nages = length(age)
        @simd for i in 1:nages
            ratio[i] = ratioPb207Pb206(age[i])
        end
        return ratio
    else
        ratio = ratioPb207Pb206(age)
    end
end

"""
    ratio2agePb206U238(ratio)

Compute the ²⁰⁶Pb/²³⁸U age for a given ratio.

Works for vectors and other arrays if column is specified.

# Examples
```julia-repl
julia> ratio2agePb206U238(0.16780392747297124)
1000.0000000000005

julia> a = [0.16780392747297124; 0.3637660130212965; 0.592611306160425];
julia> ratio2agePb206U238(a)
 1000.0000000000005
 2000.0000000000002
 2999.9999999999995
```
"""
function ratio2agePb206U238(ratio)
    if isa(ratio, Array)
        age = zeros(size(ratio))
        nratios = length(ratio)
        @simd for i in 1:nratios
            age[i] = agePb206U238(ratio[i])
        end
        return age
    else
        age = agePb206U238(ratio)
    end
end


"""
    ratio2agePb207U235(ratio)

Compute the ²⁰⁷Pb/²³⁵U age for a given ratio.
Works for vectors and other arrays if column is specified.

# Examples
```julia-repl
julia> ratio2agePb207U235(1.6774102427622641)
1000.0

julia> a = [1.6774102427622641; 6.168525608048285; 18.193083888492072];
julia> ratio2agePb207U235(a)
 1000.0
 2000.0
 3000.0
```
"""
function ratio2agePb207U235(ratio)
    if isa(ratio, Array)
        age = zeros(size(ratio))
        nratios = length(ratio)
        @simd for i in 1:nratios
            age[i] = agePb207U235(ratio[i])
        end
        return age
    else
        age = agePb207U235(ratio)
    end
end

"""
    ratio2agePb207Pb206(ratio)

Compute the ²⁰⁷Pb/²⁰⁶Pb age for a given ratio.

Works for vectors and other arrays if column is specified.
Uses the Newton-Raphson iterative method to solve for age with a tolerance of 1e-12 of the input ratio(s).

# Examples
```julia-repl
julia> ratio2agePb207Pb206(0.07253226274578425)
1000.0000000001221

julia> a = [0.07253226274578425; 0.1230419792914581; 0.22275653097919715];
julia> ratio2agePb207Pb206(a)
 1000.0000000001221
 2000.0000000046768
 3000.00000000223
```
"""
function ratio2agePb207Pb206(ratio)
    if isa(ratio, Array)
        age = zeros(size(ratio))
        nratios = length(ratio)
        @threads for i in 1:nratios
            age[i] = agePb207Pb206(ratio[i])
        end
        return age
    else
        age = agePb207Pb206(ratio)
    end
end

"""
    calcaitchisonTW(rU238Pb206, rPb207Pb206, aPb206U238, aPb207Pb206)

Computes the Aitchison distance (in Tera-Wasserburg space) for specified parameters.

Used as a measure of concordance in Tera-Wasserburg space. Unlike IsoplotR, does not multiply the value by 100.

For the Wetherill space variant see: calcaitchisonW.

# Examples
```julia-repl
julia> calcaitchisonTW(5.959336083841503, 0.07434716497717946, 1000, 1050)
0.02471346503520116
```
"""
function calcaitchisonTW(rU238Pb206, rPb207Pb206, aPb206U238, aPb207Pb206)
    if isa(rU238Pb206, Array)
        aitchdist = zeros(size(rU238Pb206))
        ncount = length(rU238Pb206)
        @simd for i in 1:ncount
            aitchdist[i] = aitchisonTW(rU238Pb206[i], rPb207Pb206[i], aPb206U238[i], aPb207Pb206[i])
        end
        return aitchdist
    else
        aitchisonTW(rU238Pb206, rPb207Pb206, aPb206U238, aPb207Pb206)
    end
end

#Base functions
function ratioPb207U235(age)
    exp(λU235 * age) - 1
end

function ratioPb206U238(age)
    exp(λU238 * age) - 1
end

function ratioPb207Pb206(age)
    ((exp(λU235 * age) - 1) / (exp(λU238 * age) - 1)) * inv(U238U235)
end

function ratioprimePb207Pb206(age) #  d/dt(((e^(λU235 * age) - 1)/(e^(λU238 * age) - 1))/137.818)
    inv(U238U235) * ((λU235) * exp(λU235 * age) * (exp(λU238 * age) - 1) - λU238 *
    exp(λU238 * age) * (exp(λU235 * age) - 1)) / (exp(λU238 * age) - 1)^2
end

function agePb206U238(ratio)
    if ratio > 0
        log(ratio + 1) / λU238
    else
        throw(ArgumentError("A negative or zero value ratio is not possible for geochemical (i.e., compositional) data.
        Please check your data."))
    end
end

function agePb207U235(ratio)
    if ratio > 0
        log(ratio + 1) / λU235
    else
        throw(ArgumentError("A negative or zero value ratio is not possible for geochemical (i.e., compositional) data.
        Please check your data."))
    end
end

function agePb207Pb206(ratio)
    if ratio <= 0
        throw(ArgumentError("A negative or zero value ratio is not possible for geochemical (i.e., compositional) data.
        Please check your data."))
    else
        t_guess = (log(ratio) + 3.21121) / 0.000586671
        ε = eps()
        n_iterations = 1
        while abs(ratio - ratioPb207Pb206(t_guess)) > eps() && n_iterations < 1e6
            t_iteration = t_guess - (ratioPb207Pb206(t_guess) - ratio) / ratioprimePb207Pb206(t_guess)
            t_guess = t_iteration
            n_iterations += 1
        end
        return if n_iterations < 1e6
            t_guess
        else
          warn(("Convergence at machine tolerance ($(eval(eps()))) not reached within 1e6 iterations!"))
        end
    end
end

function ageconcordia(rPb206U238, rPb207U235, rPb207Pb206)
    if rPb206U238 > 0.0 && rPb207U235 > 0.0 && rPb207Pb206 > 0.0
        # concordia age calculation goes here
    else
        throw(ArgumentError("A negative or zero value ratio is not possible for geochemical (i.e., compositional) data.
        Please check your data."))
    end
end

function aitchisonTW(rU238Pb206, rPb207Pb206, aPb206U238, aPb207Pb206)
    if rU238Pb206 > 0.0 && rPb207Pb206 > 0.0
        signPb = sign(aPb207Pb206)
        signU = sign(aPb206U238)
        (signPb * signU) * ((log(rU238Pb206) - log(exp(λU238 * abs(aPb207Pb206)) - 1)) * sin(atan((log(rPb207Pb206) -
        log(inv(U238U235) * (exp(λU235 * abs(aPb206U238)) - 1) / (exp(λU238 * abs(aPb206U238)) - 1))) /
        (log(rU238Pb206) - log(exp(λU238 * abs(aPb207Pb206)) - 1)))))
    else
        throw(ArgumentError("A negative or zero value ratio is not possible for geochemical (i.e., compositional) data.
        Please check your data."))
    end
end

function aitchisonW(rPb206U238, rPb207U235, aPb206U238, aPb207U235)
    if rPb206U238 > 0.0 && rPb207U235 > 0.0
        signU238a = sign(aPb206U238)
        signU235a = sign(aPb207U235)
        # Aitchison distance calculation (Wetherill space) goes here
    else
        throw(ArgumentError("A negative or zero value ratio is not possible for geochemical (i.e., compositional) data.
        Please check your data."))
    end
end
