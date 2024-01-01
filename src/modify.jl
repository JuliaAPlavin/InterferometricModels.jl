using AccessorsExtra: ConstrainedLens, FixArgs
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
             co::ConstrainedLens{<:Union{typeof(intensity_peak), Base.Fix2{typeof(Tb_peak)}, FixArgs{typeof(visibility),<:Tuple{typeof(abs),Vararg{Any}}}}, PropertyLens{:flux}},
             val)
    ratio = convert(Real, val / co(c))
    set(c, co.mo, co.mo(c) * ratio)
end

function set(c::CircularGaussian,
             co::ConstrainedLens{<:FixArgs{typeof(visibility),<:Tuple{typeof(abs),Vararg{Any}}}, PropertyLens{:σ}},
             val)
    uv = co.o.args[3]
    T = typeof(c.σ)
    val > flux(c) && error("Cannot make visibility higher than total flux")
    newval = √( log(val / c.flux) / (-2π^2 * dot(uv, uv)) )
    set(c, co.mo, convert(T, newval))
end
