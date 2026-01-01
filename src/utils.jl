σ_to_fwhm(σ) = oftype(float(σ), √(8 * log(2)) * σ)
fwhm_to_σ(fwhm) = oftype(float(fwhm), fwhm / √(8 * log(2)))
Accessors.inverse(::typeof(σ_to_fwhm)) = fwhm_to_σ
Accessors.inverse(::typeof(fwhm_to_σ)) = σ_to_fwhm


# unitful: use units
intensity_to_Tb(intensity, ν) = intensity * u"c"^2 / (2 * u"k" * ν^2) |> u"K"
Tb_to_intensity(Tb, ν) = Tb * 2 * u"k" * ν^2 / u"c"^2
# plain number: assume janskys/mas^2
intensity_to_Tb(intensity::Real, ν) = ustrip(u"K", intensity_to_Tb(intensity*u"Jy/(1e-3*arcsecond)^2", ν))
Tb_to_intensity(Tb::Real, ν) = Tb_to_intensity(Tb*u"K", ν) / u"Jy/(1e-3*arcsecond)^2" |> upreferred |> ustrip

Accessors.inverse(f::Base.Fix2{typeof(intensity_to_Tb)}) = Base.Fix2(Tb_to_intensity, f.x)
Accessors.inverse(f::Base.Fix2{typeof(Tb_to_intensity)}) = Base.Fix2(intensity_to_Tb, f.x)


const UVType = StaticVector{2}
const XYType = StaticVector{2}

position_angle(x::StaticVector{2}) = atan(first(x), last(x))

struct UVW{T} <: FieldVector{3, T}
    u::T
    v::T
    w::T
end

struct UV{T} <: FieldVector{2, T}
    u::T
    v::T
end

UV(uvw::UVW) = UV(uvw.u, uvw.v)

StaticArrays.similar_type(::Type{<:UV}, ::Type{T}, s::Size{(2,)}) where {T} = UV{T}
StaticArrays.similar_type(::Type{<:UVW}, ::Type{T}, s::Size{(3,)}) where {T} = UVW{T}


const MIN_EXP_ARG = log(1e-30)  # if flux < 10^-30 of component peak, return zero

# difmap uses lookup table for exp(), see mapres.c file in difmap source code
# here, we define exp_for_gaussintensity() that can be swapped to use either:
# - basic exp (correct, but slightly differs from difmap)
# - difmap-like exp with lookup table (when matching difmap results is important)
exp_for_gaussintensity(x) = exp_for_gaussintensity_basic(x)
exp_for_gaussintensity_basic(x) = x > MIN_EXP_ARG ? exp(x) : zero(exp(x))
function exp_for_gaussintensity_difmap(x)
    @assert x ≤ 0
    ETSIZ = 1024
    nsigma = 4.5
    expconv = ETSIZ/(0.5*nsigma^2)
    x1 = -x * expconv
    x1 > ETSIZ && return zero(exp(x))
    x2 = -floor(x1) / expconv
    return exp(x2)
end
