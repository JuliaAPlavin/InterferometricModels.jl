const UVType = StaticVector{2}
const XYType = StaticVector{2}


abstract type ModelComponent end

Base.broadcastable(c::ModelComponent) = Ref(c)

position_angle((c_from, c_to)::Pair{<:ModelComponent, <:ModelComponent}) = atan( (coords(c_to) - coords(c_from))... )
separation(c_from::ModelComponent, c_to::ModelComponent) = hypot( (coords(c_from) - coords(c_to))... )

@accessor flux(c::ModelComponent) = c.flux
@accessor coords(c::ModelComponent) = c.coords
Tb_peak(c::ModelComponent, ν) = intensity_to_Tb(intensity_peak(c), ν)


Base.@kwdef struct Point{TF,TC} <: ModelComponent
    flux::TF
    coords::SVector{2, TC}
end

fwhm_max(c::Point) = fwhm_average(c)
fwhm_min(c::Point) = fwhm_average(c)
fwhm_average(c::Point{TF,TC}) where {TF,TC} = zero(TC)
effective_area(c::Point) = fwhm_average(c)^2  # also zero, but need proper units


Base.@kwdef struct CircularGaussian{TF,TS,TC} <: ModelComponent
    flux::TF
    σ::TS
    coords::SVector{2, TC}
end

@accessor fwhm_max(c::CircularGaussian) = fwhm_average(c)
@accessor fwhm_min(c::CircularGaussian) = fwhm_average(c)
@accessor fwhm_average(c::CircularGaussian) = σ_to_fwhm(c.σ)
@accessor effective_area(c::CircularGaussian) = 2π * c.σ^2


Base.@kwdef struct EllipticGaussian{TF,TS,TC,T} <: ModelComponent
    flux::TF
    σ_major::TS
    ratio_minor_major::T
    pa_major::T
    coords::SVector{2, TC}

    function EllipticGaussian(flux, σ_major, ratio_minor_major, pa_major, coords)
        if ratio_minor_major > 1
            σ_major = (σ_major * ratio_minor_major)::typeof(σ_major)
            ratio_minor_major = inv(ratio_minor_major)::typeof(ratio_minor_major)
            pa_major = (pa_major + π/2)::typeof(pa_major)
        end
        new{typeof(flux),typeof(σ_major),eltype(coords),typeof(ratio_minor_major)}(flux, σ_major, ratio_minor_major, pa_major, coords)
    end
end

EllipticGaussian(c::EllipticGaussian) = c
EllipticGaussian(c::CircularGaussian) = EllipticGaussian(c.flux, c.σ, 1., 0., c.coords)

fwhm_max(c::EllipticGaussian) = σ_to_fwhm(c.σ_major)
fwhm_min(c::EllipticGaussian) = σ_to_fwhm(c.σ_major * c.ratio_minor_major)
fwhm_average(c::EllipticGaussian) = σ_to_fwhm(c.σ_major * √(c.ratio_minor_major))  # geometric average
effective_area(c::EllipticGaussian) = 2π * c.σ_major^2 * c.ratio_minor_major
@accessor position_angle(c::EllipticGaussian) = c.pa_major
Accessors.set(c::EllipticGaussian, ::typeof(fwhm_max), val) = setproperties(c, (σ_major = fwhm_to_σ(val), ratio_minor_major=c.σ_major/fwhm_to_σ(val) * c.ratio_minor_major))
Accessors.set(c::EllipticGaussian, ::typeof(fwhm_min), val) = @set c.ratio_minor_major = fwhm_to_σ(val) / c.σ_major


Base.@kwdef struct EllipticGaussianCovmat{TF,TC,TM} <: ModelComponent
    flux::TF
    covmat::TM
    coords::SVector{2, TC}
end

fwhm_max(c::EllipticGaussianCovmat) = √(eigen(c.covmat).values[2]) |> σ_to_fwhm
fwhm_min(c::EllipticGaussianCovmat) = √(eigen(c.covmat).values[1]) |> σ_to_fwhm
fwhm_average(c::EllipticGaussianCovmat) = (eigen(c.covmat).values |> prod)^0.25 |> σ_to_fwhm

effective_area(c::EllipticGaussianCovmat) = 2π * sqrt(det(c.covmat))

EllipticGaussianCovmat(c::EllipticGaussianCovmat) = c
EllipticGaussianCovmat(c::CircularGaussian) = EllipticGaussianCovmat(; c.flux, covmat=Diagonal(SVector(1, 1) .* c.σ^2), c.coords)
EllipticGaussianCovmat(c::EllipticGaussian) = let
    si, co = sincos(c.pa_major)
    σs = c.σ_major .* SVector(c.ratio_minor_major, 1)
    trmat = @SMatrix([co si; -si co]) * Diagonal(σs)
    covmat = trmat * trmat'
    EllipticGaussianCovmat(; c.flux, covmat, c.coords)
end

function EllipticGaussian(c::EllipticGaussianCovmat)
    E = eigen(c.covmat)
    vec = E.vectors[:, 2]
    EllipticGaussian(;
        c.flux, c.coords,
        σ_major=√(E.values[2]),
        ratio_minor_major=√(E.values[1] / E.values[2]),
        pa_major=atan(vec[1], vec[2]),
    )
end


Base.@kwdef struct MultiComponentModel{TUP}
    components::TUP
end

Base.broadcastable(c::MultiComponentModel) = Ref(c)
@accessor components(m::MultiComponentModel) = m.components

Base.:(==)(a::MultiComponentModel, b::MultiComponentModel) = components(a) == components(b)

flux(m::MultiComponentModel) = sum(flux, components(m))

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
    invcovmat = inv(c.covmat)
    (xy::XYType) -> let
        Δxy = xy - c.coords
        peak * exp(-1/2 * dot(Δxy, invcovmat, Δxy))
    end
end
function intensity(m::MultiComponentModel)
    bycomp = components(m) .|> intensity
    (xy::XYType) -> sum(f -> f(xy), bycomp)
end


# default definitions: work for all symmetric components
@inline visibility(::typeof(angle), c::ModelComponent, uv::UVType) = 2π * dot(uv, coords(c))
@inline visibility_envelope(::typeof(angle), c::ModelComponent, uvdist::Real) = 0 ± 2π * min(norm(coords(c)) * uvdist, 0.5)
@inline visibility_envelope(f::ComposedFunction{<:Any, typeof(angle)}, c, uvdist::Real) =
    @modify(visibility_envelope(angle, c, uvdist) |> endpoints |> _[∗]) do x
        # XXX: assuming f is monotonic wrt angle
        f.outer(x)
    end

visibility(::typeof(abs), c::Point, uv::UVType) = c.flux
visibility(::typeof(abs), c::CircularGaussian, uv::UVType) = c.flux * exp(-2π^2 * c.σ^2 * dot(uv, uv))
visibility(::typeof(abs), c::EllipticGaussian, uv::UVType) = visibility(abs, EllipticGaussianCovmat(c), uv)
visibility(::typeof(abs), c::EllipticGaussianCovmat, uv::UVType) = c.flux * exp(-2π^2 * dot(uv, c.covmat, uv))

visibility(c::Union{ModelComponent,MultiComponentModel}, uv::UVType) = visibility(c)(uv)
visibility(c::Point) = @inline (uv::UVType) -> flux(c) * cis(visibility(angle, c, uv))
function visibility(c::CircularGaussian)
    mul = -2π^2 * c.σ^2
    @inline (uv::UVType) -> flux(c) * exp(mul * dot(uv, uv)) * cis(visibility(angle, c, uv))
end
visibility(c::EllipticGaussian) = visibility(EllipticGaussianCovmat(c))
visibility(c::EllipticGaussianCovmat) = @inline uv::UVType -> flux(c) * exp(-2π^2 * dot(uv, c.covmat * uv)) * cis(visibility(angle, c, uv))
function visibility(m::MultiComponentModel)
    bycomp = components(m) .|> visibility
    @inline (xy::UVType) -> sum(f -> f(xy), bycomp)
end

visibility_envelope(::typeof(abs), c::Point, uvdist::Real) = c.flux ± eps(float(c.flux))
visibility_envelope(::typeof(abs), c::CircularGaussian, uvdist::Real) = c.flux * exp(-2π^2 * c.σ^2 * uvdist^2) ± eps(float(c.flux))
visibility_envelope(::typeof(abs), c::EllipticGaussian, uvdist::Real) = (c.flux * exp(-2π^2 * c.σ_major^2 * uvdist^2)) .. (c.flux * exp(-2π^2 * (c.σ_major * c.ratio_minor_major)^2 * uvdist^2))

visibility_envelope(f::typeof(abs), m::MultiComponentModel{<:Tuple{Any}}, uvdist::Real) = visibility_envelope(f, only(components(m)), uvdist)
visibility_envelope(f::typeof(angle), m::MultiComponentModel{<:Tuple{Any}}, uvdist::Real) = visibility_envelope(f, only(components(m)), uvdist)

visibility(f::Function, args...; kwargs...) = f(visibility(args...; kwargs...))
Broadcast.broadcasted(::typeof(visibility), f::Function, args...; kwargs...) = f.(visibility.(args...; kwargs...))


Unitful.ustrip(x::ModelComponent) = @modify(x -> ustrip.(x), x |> Properties())
Unitful.ustrip(x::MultiComponentModel) = @modify(ustrip, x.components |> Elements())
# piracy, but...
Unitful.ustrip(x::AbstractInterval) = @modify(ustrip, x |> Properties())
Unitful.ustrip(u::Unitful.Units, x::AbstractInterval) = @modify(f -> ustrip(u, f), x |> Properties())


Base.isapprox(a::ModelComponent, b::ModelComponent; kwargs...) =
    propertynames(a) == propertynames(b) && all(k -> isapprox(getproperty(a, k), getproperty(b, k); kwargs...), propertynames(a))
Base.isapprox(a::MultiComponentModel, b::MultiComponentModel; kwargs...) =
    length(components(a)) == length(components(b)) && all(((x, y),) -> isapprox(x, y; kwargs...), zip(components(a), components(b)))


import AccessorsExtra

AccessorsExtra.construct(::Type{Point}, (_,flux)::Pair{typeof(flux)}, (_,coords)::Pair{typeof(coords)}) =
    Point(; flux, coords=SVector{2}(coords))
AccessorsExtra.construct(::Type{CircularGaussian}, (_,flux)::Pair{typeof(flux)}, (_,σ)::Pair{PropertyLens{:σ}}, (_,coords)::Pair{typeof(coords)}) =
    CircularGaussian(; flux, σ=σ, coords=SVector{2}(coords))
AccessorsExtra.construct(::Type{CircularGaussian}, (_,flux)::Pair{typeof(flux)}, (_,fwhm)::Pair{typeof(fwhm_average)}, (_,coords)::Pair{typeof(coords)}) =
    CircularGaussian(; flux, σ=fwhm_to_σ(fwhm), coords=SVector{2}(coords))
AccessorsExtra.construct(::Type{EllipticGaussian}, (_,flux)::Pair{typeof(flux)}, (_,fwhm_max)::Pair{typeof(fwhm_max)}, (_,ratio)::Pair{PropertyLens{:ratio_minor_major}}, (_,pa)::Pair{typeof(position_angle)}, (_,coords)::Pair{typeof(coords)}) =
    EllipticGaussian(; flux, σ_major=fwhm_to_σ(fwhm_max), ratio_minor_major=ratio, pa_major=pa, coords=SVector{2}(coords))
AccessorsExtra.construct(::Type{EllipticGaussian}, (_,flux)::Pair{typeof(flux)}, (_,fwhm_max)::Pair{typeof(fwhm_max)}, (_,fwhm_min)::Pair{typeof(fwhm_min)}, (_,pa)::Pair{typeof(position_angle)}, (_,coords)::Pair{typeof(coords)}) =
    EllipticGaussian(; flux, σ_major=fwhm_to_σ(fwhm_max), ratio_minor_major=fwhm_min/fwhm_max, pa_major=pa, coords=SVector{2}(coords))
