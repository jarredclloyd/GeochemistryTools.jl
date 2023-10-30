using LinearAlgebra
function covmat_rho(sx, sy, rxy; se_level_in = 2)
    sxy = @. rxy * (sx /se_level_in) * (sy / se_level_in)
    A =@. [(sx /se_level_in)^2 sxy; sxy (sy / se_level_in)^2]
    return A
end

function getellipsepoints(centre, xradius, yradius, θ)
    t = range(0, 2 * pi; length = 100)
    ellipse_x_r = @. xradius * cos(t)
    ellipse_y_r = @. yradius * sin(t)
    R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
    r_ellipse = [ellipse_x_r ellipse_y_r] * R
    x = @. centre[1] + r_ellipse[:, 1]
    y = @. centre[2] + r_ellipse[:, 2]
    return zip(x,y)
end

function getellipsepoints(centre, Σ, confidence = 0.95)
    quant = quantile(Chisq(2), confidence) |> sqrt
    smallesttegv = eigmin(Σ)
    largestegv = eigmax(Σ)
    idxmax = findmax(eigvals(Σ))[2]

    rx = quant * sqrt(largestegv)
    ry = quant * sqrt(smallesttegv)

    eigvecmax = eigvecs(Σ)[:, idxmax]
    θ = atan(eigvecmax[2] / eigvecmax[1])
    if θ < 0
        θ += 2 * π
    end

    return getellipsepoints(centre, rx, ry, θ)
end

function errorellipse(cx, cy, sx, sy, rxy; se_level_in = 2, confidence = 0.95)
    return getellipsepoints([cx, cy], covmat_rho(sx, sy, rxy; se_level_in = se_level_in), confidence)
end
