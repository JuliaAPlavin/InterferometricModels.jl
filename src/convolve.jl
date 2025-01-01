function convolve(c::Point, b::CircularGaussian)
    @assert coords(b) == SVector(0, 0)
    return CircularGaussian(
        flux=flux(c) * flux(b),
        σ=b.σ,
        coords=coords(c),
    )
end

function convolve(c::Point, b::EllipticGaussian)
    @assert coords(b) == SVector(0, 0)
    return EllipticGaussian(;
        flux=flux(c) * flux(b),
        b.σ_major,
        b.ratio_minor_major,
        b.pa_major,
        coords=coords(c),
    )
end

function convolve(c::CircularGaussian, b::CircularGaussian)
    @assert coords(b) == SVector(0, 0)
    return CircularGaussian(
        flux=flux(c) * flux(b),
        σ=hypot(c.σ, b.σ),
        coords=coords(c),
    )
end

let TS = [CircularGaussian, EllipticGaussian, EllipticGaussianCovmat]
    @eval convolve(c::Union{$TS...}, b::Union{$TS...}) = convolve(EllipticGaussianCovmat(c), EllipticGaussianCovmat(b))
end
function convolve(c::EllipticGaussianCovmat, b::EllipticGaussianCovmat)
    @assert coords(b) == SVector(0, 0)
    return EllipticGaussianCovmat(
        flux=flux(c) * flux(b),
        covmat=c.covmat + b.covmat,
        coords=coords(c),
    )
end

convolve(m::MultiComponentModel, b::ModelComponent) = MultiComponentModel(convolve.(components(m), b))
