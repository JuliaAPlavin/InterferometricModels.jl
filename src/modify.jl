# type piracy:
# doesn't change existing behavior, just extends to more than two arguments
(f::Base.Fix1)(args...) = f.f(f.x, args...)
(f::Base.Fix2)(y, args...) = f.f(y, f.x, args...)


Accessors.set(c::ModelComponent, ::typeof(flux), val) = @set c.flux = val
Accessors.set(c::ModelComponent, ::typeof(coords), val) = @set c.coords = val
Accessors.set(c::CircularGaussian, ::Union{typeof(fwhm_max), typeof(fwhm_min), typeof(fwhm_average)}, val) = @set c.σ = fwhm_to_σ(val)
Accessors.set(c::EllipticGaussian, ::typeof(fwhm_max), val) = setproperties(c, (σ_major = fwhm_to_σ(val), ratio_minor_major=c.σ_major/fwhm_to_σ(val) * c.ratio_minor_major))
Accessors.set(c::EllipticGaussian, ::typeof(fwhm_min), val) = @set c.ratio_minor_major = fwhm_to_σ(val) / c.σ_major
Accessors.set(c::MultiComponentModel, ::typeof(components), val) = @set c.components = val

"""    set_so_that(x, optic, func => target_val)

Set `optic` in `x` so that `func(new_x) == target_val`.
"""
set_so_that(x, o, (func, val)::Pair) = set_so_that(x, o, func, val)

function set_so_that(
        c::CircularGaussian, o::Accessors.PropertyLens{:σ},
        f::Union{typeof(fwhm_max), typeof(fwhm_min), typeof(fwhm_average)},
        val)
    ratio = convert(Real, val / f(c))
    set(c, o, o(c) * ratio)
end

function set_so_that(
        c::CircularGaussian, o::Accessors.PropertyLens{:σ},
        f::Union{typeof(effective_area)},
        val)
    ratio = convert(Real, val / f(c))
    set(c, o, o(c) * √ratio)
end

function set_so_that(
        c::CircularGaussian, o::Accessors.PropertyLens{:σ},
        f::Union{typeof(intensity_peak), Base.Fix2{typeof(Tb_peak)}},
        val)
    ratio = convert(Real, val / f(c))
    set(c, o, o(c) / √(ratio))
end

function set_so_that(
        c::CircularGaussian, o::Accessors.PropertyLens{:flux},
        f::Union{typeof(intensity_peak), Base.Fix2{typeof(Tb_peak)}, Base.Fix2{Base.Fix1{typeof(visibility), typeof(abs)}}},
        val)
    ratio = convert(Real, val / f(c))
    set(c, o, o(c) * ratio)
end

function set_so_that(
        c::CircularGaussian, o::Accessors.PropertyLens{:σ},
        f::Base.Fix2{Base.Fix1{typeof(visibility), typeof(abs)}},
        val)
    uv = f.x
    T = typeof(c.σ)
    val > flux(c) && error("Cannot make visibility higher than total flux")
    newval = √( log(val / c.flux) / (-2π^2 * dot(uv, uv)) )
    @set c.σ = convert(T, newval)
end
