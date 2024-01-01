# type piracy:
# doesn't change existing behavior, just extends to more than two arguments
(f::Base.Fix1)(args...) = f.f(f.x, args...)
(f::Base.Fix2)(y, args...) = f.f(y, f.x, args...)

using AccessorsExtra: ConstrainedLens
import Accessors: set

function set(c::CircularGaussian,
             co::ConstrainedLens{<:Union{typeof(fwhm_max), typeof(fwhm_min), typeof(fwhm_average)}, PropertyLens{:σ}},
             val)
    ratio = convert(Real, val / co(c))
    set(c, co.mo, co.mo(c) * ratio)
end

function set(c::CircularGaussian,
             co::ConstrainedLens{typeof(effective_area), PropertyLens{:σ}},
             val)
    ratio = convert(Real, val / co(c))
    set(c, co.mo, co.mo(c) * √ratio)
end

function set(c::CircularGaussian,
             co::ConstrainedLens{<:Union{typeof(intensity_peak), Base.Fix2{typeof(Tb_peak)}}, PropertyLens{:σ}},
             val)
    ratio = convert(Real, val / co(c))
    set(c, co.mo, co.mo(c) / √ratio)
end

function set(c::CircularGaussian,
             co::ConstrainedLens{<:Union{typeof(intensity_peak), Base.Fix2{typeof(Tb_peak)}, Base.Fix2{Base.Fix1{typeof(visibility), typeof(abs)}}}, PropertyLens{:flux}},
             val)
    ratio = convert(Real, val / co(c))
    set(c, co.mo, co.mo(c) * ratio)
end

function set(c::CircularGaussian,
             co::ConstrainedLens{<:Base.Fix2{Base.Fix1{typeof(visibility), typeof(abs)}}, PropertyLens{:σ}},
             val)
    uv = co.o.x
    T = typeof(c.σ)
    val > flux(c) && error("Cannot make visibility higher than total flux")
    newval = √( log(val / c.flux) / (-2π^2 * dot(uv, uv)) )
    set(c, co.mo, convert(T, newval))
end
