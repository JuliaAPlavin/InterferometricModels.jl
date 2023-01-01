const UVType = StaticVector{2}
const XYType = StaticVector{2}


abstract type ModelComponent end

Base.broadcastable(c::ModelComponent) = Ref(c)

position_angle((c_from, c_to)::Pair{<:ModelComponent, <:ModelComponent}) = atan( (coords(c_to) - coords(c_from))... )
separation(c_from::ModelComponent, c_to::ModelComponent) = hypot( (coords(c_from) - coords(c_to))... )

flux(c::ModelComponent) = c.flux
coords(c::ModelComponent) = c.coords
Tb_peak(c::ModelComponent, ν) = intensity_to_Tb(intensity_peak(c), ν)


Base.@kwdef struct Point{TF,TC} <: ModelComponent
    flux::TF
    coords::SVector{2, TC}
end

fwhm_max(c::Point) = fwhm_average(c)
fwhm_min(c::Point) = fwhm_average(c)
fwhm_average(c::Point{TF,TC}) where {TF,TC} = zero(TC)
effective_area(c::Point) = fwhm_average(c)^2  # also zero, but need proper units


Base.@kwdef struct CircularGaussian{TF,TC} <: ModelComponent
    flux::TF
    σ::TC
    coords::SVector{2, TC}
end

fwhm_max(c::CircularGaussian) = fwhm_average(c)
fwhm_min(c::CircularGaussian) = fwhm_average(c)
fwhm_average(c::CircularGaussian) = σ_to_fwhm(c.σ)
effective_area(c::CircularGaussian) = 2π * c.σ^2


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
position_angle(c::EllipticGaussian) = c.pa_major


Base.@kwdef struct EllipticGaussianCovmat{TF,TC,TM} <: ModelComponent
    flux::TF
    covmat::TM
    coords::SVector{2, TC}
end

effective_area(c::EllipticGaussianCovmat) = 2π * sqrt(1/det(c.covmat))

EllipticGaussianCovmat(c::CircularGaussian) = EllipticGaussianCovmat(; c.flux, covmat=Diagonal(SVector(1, 1) ./ c.σ^2), c.coords)
EllipticGaussianCovmat(c::EllipticGaussian) = let
    si, co = sincos(position_angle(c))
    σs_inv = 1 ./ (c.σ_major .* SVector(c.ratio_minor_major, 1))
    trmat = Diagonal(σs_inv) * @SMatrix([co -si; si co])
    covmat = trmat' * trmat
    EllipticGaussianCovmat(; c.flux, covmat, c.coords)
end


Base.@kwdef struct MultiComponentModel{TUP}
    components::TUP
end

Base.broadcastable(c::MultiComponentModel) = Ref(c)
components(m::MultiComponentModel) = m.components


intensity_peak(c::Point) = flux(c) / effective_area(c)
intensity_peak(c::CircularGaussian) = flux(c) / effective_area(c)
intensity_peak(c::EllipticGaussian) = flux(c) / effective_area(c)
intensity_peak(c::EllipticGaussianCovmat) = flux(c) / effective_area(c)

function intensity(c::CircularGaussian)
    peak = intensity_peak(c)
    mul = -1 / (2*c.σ^2)
    (xy::XYType) -> peak * exp(mul * dot(xy - c.coords, xy - c.coords))
end
intensity(c::EllipticGaussian) = intensity(EllipticGaussianCovmat(c))
function intensity(c::EllipticGaussianCovmat)
    peak = intensity_peak(c)
    (xy::XYType) -> let
        Δxy = xy - c.coords
        peak * exp(-1/2 * dot(Δxy, c.covmat, Δxy))
    end
end
function intensity(m::MultiComponentModel)
    bycomp = components(m) .|> intensity
    (xy::XYType) -> sum(i -> i(xy), bycomp)
end


# default definitions: work for all symmetric components
@inline visibility(::typeof(angle), c::ModelComponent, uv::UVType) = 2π * dot(uv, coords(c))
@inline visibility_envelope(::typeof(angle), c::ModelComponent, uvdist::Real) = 0 ± 2π * min(norm(coords(c)) * uvdist, 0.5)

visibility(::typeof(abs), c::Point, uv::UVType) = c.flux
visibility(::typeof(abs), c::CircularGaussian, uv::UVType) = c.flux * exp(-2π^2 * c.σ^2 * dot(uv, uv))
visibility(::typeof(abs), c::EllipticGaussian, uv::UVType) = visibility(abs, EllipticGaussianCovmat(c), uv)
visibility(::typeof(abs), c::EllipticGaussianCovmat, uv::UVType) = c.flux * exp(-2π^2 * dot(uv, inv(c.covmat), uv))

visibility(c::Point, uv::UVType) = c.flux * exp(im * visibility(angle, c, uv))
# with fwhm: c.flux * exp(-pi^2 / log(16) * c.fwhm^2 * dot(uv, uv) - 2π*im * dot(uv, c.coords))
visibility(c::CircularGaussian, uv::UVType) = c.flux * exp(-2π^2 * c.σ^2 * dot(uv, uv) + im * visibility(angle, c, uv))
visibility(c::EllipticGaussian, uv::UVType) = visibility(EllipticGaussianCovmat(c), uv)
visibility(c::EllipticGaussianCovmat, uv::UVType) = c.flux * exp(-2π^2 * dot(uv, inv(c.covmat), uv) + im * visibility(angle, c, uv))
visibility(m::MultiComponentModel, uv::UVType) = sum(c -> visibility(c, uv), components(m))

visibility_envelope(::typeof(abs), c::Point, uvdist::Real) = c.flux ± eps(c.flux)
visibility_envelope(::typeof(abs), c::CircularGaussian, uvdist::Real) = c.flux * exp(-2π^2 * c.σ^2 * uvdist^2) ± eps(c.flux)
visibility_envelope(::typeof(abs), c::EllipticGaussian, uvdist::Real) = (c.flux * exp(-2π^2 * c.σ_major^2 * uvdist^2)) .. (c.flux * exp(-2π^2 * (c.σ_major * c.ratio_minor_major)^2 * uvdist^2))


Unitful.ustrip(x::ModelComponent) = @modify(ustrip, x |> Properties())
Unitful.ustrip(x::MultiComponentModel) = @modify(ustrip, x.components |> Elements())
# piracy, but...
Unitful.ustrip(x::AbstractInterval) = @modify(ustrip, x |> Properties())
Unitful.ustrip(u::Unitful.Units, x::AbstractInterval) = @modify(f -> ustrip(u, f), x |> Properties())


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
        covmat=inv(inv(c.covmat) + inv(b.covmat)),
        coords=coords(c),
    )
end

convolve(m::MultiComponentModel, b::ModelComponent) = MultiComponentModel(convolve.(components(m), b))
