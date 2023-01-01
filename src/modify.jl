# type piracy:
# doesn't change existing behavior, just extends to more than two arguments
(f::Base.Fix1)(args...) = f.f(f.x, args...)
(f::Base.Fix2)(y, args...) = f.f(y, f.x, args...)

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
    newval = √( log(val / c.flux) / (-2π^2 * dot(uv, uv)) )
    @set c.σ = convert(T, newval)
end
