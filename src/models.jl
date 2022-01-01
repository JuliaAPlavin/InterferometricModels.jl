const UVType = StaticVector{2}
const XYType = StaticVector{2}


module Transformations

abstract type Transformation end
function transform end

struct Identity <: Transformation end
transform(tr::Identity, x, ::Val) = x

Base.@kwdef struct FuncTransformation{T1, T2, T3} <: Transformation
	flux::T1
	σ::T2
	x::T3
end
transform(tr::FuncTransformation, x, ::Val{:flux}) = tr.flux(x)
transform(tr::FuncTransformation, x, ::Val{:σ}) = tr.σ(x)
transform(tr::FuncTransformation, x, ::Val{:coords}) = (tr.x(x[1]), tr.x(x[2]))
end


abstract type ModelComponent end

position_angle((c_from, c_to)::Pair{<:ModelComponent, <:ModelComponent}) = atan( (coords(c_to) - coords(c_from))... )
distance(c_from::ModelComponent, c_to::ModelComponent) = hypot( (coords(c_from) - coords(c_to))... )

flux(c::ModelComponent) = c.flux
coords(c::ModelComponent) = c.coords
Tb_peak(c::ModelComponent, ν) = intensity_to_Tb(intensity_peak(c), ν)

# unitful: use units
intensity_to_Tb(intensity, ν) = intensity * u"c"^2 / (2 * u"k" * ν^2) |> u"K"
# plain number: assume janskys/mas^2
intensity_to_Tb(intensity::Real, ν) = ustrip(u"K", intensity_to_Tb(intensity*u"Jy/(1e-3*arcsecond)^2", ν))


Base.@kwdef struct PointSource{TF,TC} <: ModelComponent
    flux::TF
    coords::SVector{2, TC}
end

fwhm_max(c::PointSource) = fwhm_average(c)
fwhm_min(c::PointSource) = fwhm_average(c)
fwhm_average(c::PointSource{TF,TC}) where {TF,TC} = zero(TC)
effective_area(c::PointSource) = fwhm_average(c)^2  # also zero, but need proper units
intensity_peak(c::PointSource) = flux(c) / effective_area(c)

visibility_amplitude(c::PointSource, uv::UVType) = c.flux
visibility_phase(c::PointSource, uv::UVType) = 2π * dot(uv, c.coords)
function visibility(c::PointSource, uv::UVType)
    c.flux * exp(-π2im * dot(uv, c.coords))
end
visibility_envelope(c::PointSource, uvdist::Real) = c.flux ± zero(typeof(c.flux))


Base.@kwdef struct CircularGaussian{TF,TC} <: ModelComponent
    flux::TF
    σ::TC
    coords::SVector{2, TC}
end

n_params(::Type{<:CircularGaussian}) = 4
function from_vec(::Type{<:CircularGaussian}, x::AbstractVector, start=1)
    CircularGaussian(x[start], x[start+1], SVector(x[start+2], x[start+3]))
end
function transform_vec!(::Type{<:CircularGaussian}, tr::Transformations.Transformation, x::AbstractVector, start=1)
    x[start] = Transformations.transform(tr, x[start], Val(:flux))
    x[start+1] = Transformations.transform(tr, x[start+1], Val(:σ))
    x[start+2], x[start+3] = Transformations.transform(tr, SVector(x[start+2], x[start+3]), Val(:coords))
end
param_names(::Type{<:CircularGaussian}) = ["flux", "σ", "x", "y"]

fwhm_max(c::CircularGaussian) = fwhm_average(c)
fwhm_min(c::CircularGaussian) = fwhm_average(c)
fwhm_average(c::CircularGaussian) = σ_to_fwhm(c.σ)
effective_area(c::CircularGaussian) = 2π * c.σ^2

intensity_peak(c::CircularGaussian) = flux(c) / effective_area(c)
intensity(c::CircularGaussian, xy::XYType) = intensity_peak(c) * exp(-dot(xy - c.coords, xy - c.coords)/(2*c.σ^2))

visibility_amplitude(c::CircularGaussian, uv::UVType) = c.flux * exp(-2π^2 * c.σ^2 * dot(uv, uv))
visibility_phase(c::CircularGaussian, uv::UVType) = 2π * dot(uv, c.coords)
function visibility(c::CircularGaussian, uv::UVType)
    c.flux * exp(-2π^2 * c.σ^2 * dot(uv, uv) - 2π*im * dot(uv, c.coords))
    # with fwhm instead of σ:
    # c.flux * exp(-pi^2 / log(16) * c.fwhm^2 * dot(uv, uv) - 2π*im * dot(uv, c.coords))
end
visibility_envelope(c::CircularGaussian, uvdist::Real) = c.flux * exp(-2π^2 * c.σ^2 * uvdist^2) ± 0


Base.@kwdef struct EllipticGaussian{TF,TC,T} <: ModelComponent
    flux::TF
    σ_major::TC
    ratio_minor_major::T
    pa_major::T
    coords::SVector{2, TC}
end

fwhm_max(c::EllipticGaussian) = σ_to_fwhm(c.σ_major)
fwhm_min(c::EllipticGaussian) = σ_to_fwhm(c.σ_major * c.ratio_minor_major)
fwhm_average(c::EllipticGaussian) = σ_to_fwhm(c.σ_major * √(c.ratio_minor_major))  # geometric average
effective_area(c::EllipticGaussian) = 2π * c.σ_major^2 * c.ratio_minor_major
intensity_peak(c::EllipticGaussian) = flux(c) / effective_area(c)
visibility_phase(c::EllipticGaussian, uv::UVType) = 2π * dot(uv, c.coords)
visibility_envelope(c::EllipticGaussian, uvdist::Real) = (c.flux * exp(-2π^2 * c.σ_major^2 * uvdist^2)) .. (c.flux * exp(-2π^2 * (c.σ_major * c.ratio_minor_major)^2 * uvdist^2))
position_angle(c::EllipticGaussian) = c.pa_major


Base.@kwdef struct MultiComponentModel{TUP}
    components::TUP
end

components(m::MultiComponentModel) = m.components
visibility(model::MultiComponentModel, uv::UVType) = sum(c -> visibility(c, uv), components(model))

@generated n_params(::Type{MultiComponentModel{TUP}}) where {TUP <: Tuple} = sum(n_params(T) for T in fieldtypes(TUP))
@generated function from_vec(::Type{MultiComponentModel{TUP}}, x::AbstractVector, start=1) where {TUP <: Tuple}
    types = fieldtypes(TUP)
    quote
        @assert length(x) == $(n_params(MultiComponentModel{TUP}))
        MultiComponentModel(tuple(
            $([:(from_vec($T, x, start + $(sum(n_params(TT) for TT in types[1:i-1]; init=0)))) for (i, T) in enumerate(types)]...)
        ))
    end
end

@generated function transform_vec!(::Type{MultiComponentModel{TUP}}, tr::Transformations.Transformation, x::AbstractVector, start=1) where {TUP <: Tuple}
    types = fieldtypes(TUP)
    quote
        @assert length(x) == $(n_params(MultiComponentModel{TUP}))
        $([:(transform_vec!($T, tr, x, start + $(sum(n_params(TT) for TT in types[1:i-1]; init=0)))) for (i, T) in enumerate(types)]...)
        return x
    end
end
transform_vec(TMOD::Type{<:MultiComponentModel}, tr::Transformations.Transformation, x::AbstractVector) = transform_vec!(TMOD, tr, copy(x))

param_names(::Type{MultiComponentModel{TUP}}) where {TUP <: Tuple} = mapreduce(vcat, fieldtypes(TUP) |> enumerate) do (i, T)
    ["$(i)_$(n)" for n in param_names(T)]
end


function model_from_difmap(components::AbstractVector)
    comps = map(components) do c
        coords = c.radec
        flux = c.flux
        if c.major == 0
            PointSource(; flux, coords)
        elseif c.ratio == 1
            CircularGaussian(; flux, σ=fwhm_to_σ(c.major), coords)
        else
            EllipticGaussian(; flux, σ_major=fwhm_to_σ(c.major), ratio_minor_major=c.ratio, pa_major=c.phi, coords)
        end
    end
    MultiComponentModel(Tuple(comps))
end
