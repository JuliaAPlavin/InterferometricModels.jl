using AccessorsExtra: ConstrainedLens, FixArgs

set(c::CircularGaussian,
    co::ConstrainedLens{<:Union{typeof(fwhm_max), typeof(fwhm_min), typeof(fwhm_average), typeof(effective_area)}, PropertyLens{:σ}},
    val) = set(c, co.o, val)

function set(c::CircularGaussian,
             co::ConstrainedLens{<:Union{typeof(intensity_peak), Base.Fix2{typeof(Tb_peak)}}, PropertyLens{:σ}},
             val)
    c_tmp = set(c, co.mo, oneunit(co.mo(c)))
    ratio = convert(Real, val / co(c_tmp))
    set(c_tmp, co.mo, co.mo(c_tmp) / √ratio)
end

function set(c::Union{CircularGaussian,EllipticGaussian},
             co::ConstrainedLens{<:Union{typeof(intensity_peak), Base.Fix2{typeof(Tb_peak)}, FixArgs{typeof(visibility),<:Tuple{typeof(abs),Vararg{Any}}}}, PropertyLens{:flux}},
             val)
    c_tmp = set(c, co.mo, oneunit(co.mo(c)))
    ratio = convert(Real, val / co(c_tmp))
    set(c_tmp, co.mo, co.mo(c_tmp) * ratio)
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

set(c::EllipticGaussian, ::ConstrainedLens{typeof(fwhm_max),PropertyLens{:σ_major}}, v) = @set c.σ_major = InterferometricModels.fwhm_to_σ(v)
