
is_beam(c::ModelComponent) = coords(c) == SVector(0, 0) && intensity_peak(c) == 1
beam(::Type{CircularGaussian}; σ) = CircularGaussian(; flux=2π*σ^2, σ, coords=SVector(zero(σ), zero(σ)))
beam(::Type{EllipticGaussian}; σ_major, ratio_minor_major, pa_major) =
    EllipticGaussian(; flux=2π*σ_major^2*ratio_minor_major, σ_major, ratio_minor_major, pa_major, coords=SVector(zero(σ_major), zero(σ_major)))


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
        covmat=inv(c.invcovmat + b.invcovmat),
        invcovmat=c.invcovmat + b.invcovmat,
        coords=coords(c),
    )
end

convolve(m::MultiComponentModel, b::ModelComponent) = MultiComponentModel(convolve.(components(m), b))
