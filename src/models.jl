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
@batteries Point
Point(flux, coords) = Point(flux, SVector{2}(coords))

fwhm_max(c::Point) = fwhm_average(c)
fwhm_min(c::Point) = fwhm_average(c)
fwhm_average(c::Point{TF,TC}) where {TF,TC} = zero(TC)
effective_area(c::Point) = fwhm_average(c)^2  # also zero, but need proper units


Base.@kwdef struct CircularGaussian{TF,TS,TC} <: ModelComponent
    flux::TF
    σ::TS
    coords::SVector{2, TC}
end
@batteries CircularGaussian
CircularGaussian(flux, σ, coords) = CircularGaussian(flux, σ, SVector{2}(coords))
function CircularGaussian(c::Point)
    FT = promote_type(typeof(c.flux), eltype(c.coords))
    CircularGaussian(c.flux, zero(FT), c.coords)
end
function Point(c::CircularGaussian)
    iszero(c.σ) || throw(ArgumentError("Cannot convert CircularGaussian with non-zero σ ($(c.σ)) to Point"))
    Point(c.flux, c.coords)
end

@accessor fwhm_max(c::CircularGaussian) = fwhm_average(c)
@accessor fwhm_min(c::CircularGaussian) = fwhm_average(c)
@accessor fwhm_average(c::CircularGaussian) = σ_to_fwhm(c.σ)
@accessor effective_area(c::CircularGaussian) = 2π * square(c.σ)


Base.@kwdef struct EllipticGaussian{TF,TS,TC,TR,TP} <: ModelComponent
    flux::TF
    σ_major::TS
    ratio_minor_major::TR
    pa_major::TP
    coords::SVector{2, TC}

    function EllipticGaussian(flux, σ_major, ratio_minor_major, pa_major, coords)
        if ratio_minor_major isa Number && ratio_minor_major > 1
            σ_major = (σ_major * ratio_minor_major)::typeof(σ_major)
            ratio_minor_major = inv(ratio_minor_major)::typeof(ratio_minor_major)
            pa_major = (pa_major + π/oftype(pa_major, 2))::typeof(pa_major)
        end
        new{typeof(flux),typeof(σ_major),eltype(coords),typeof(ratio_minor_major),typeof(pa_major)}(flux, σ_major, ratio_minor_major, pa_major, coords)
    end
end
@batteries EllipticGaussian

function EllipticGaussian(c::Point)
	FT = promote_type(typeof(c.flux), eltype(c.coords))
	EllipticGaussian(c.flux, zero(eltype(c.coords)), one(FT), zero(FT), c.coords)
end
function EllipticGaussian(c::CircularGaussian)
	FT = promote_type(typeof(c.flux), typeof(c.σ), eltype(c.coords))
	EllipticGaussian(c.flux, c.σ, one(FT), zero(FT), c.coords)
end
function Point(c::EllipticGaussian)
    iszero(c.σ_major) || throw(ArgumentError("Cannot convert EllipticGaussian with non-zero σ_major ($(c.σ_major)) to Point"))
    Point(c.flux, c.coords)
end
function CircularGaussian(c::EllipticGaussian)
    isone(c.ratio_minor_major) || throw(ArgumentError("Cannot convert EllipticGaussian with non-unity ratio_minor_major ($(c.ratio_minor_major)) to CircularGaussian"))
    CircularGaussian(c.flux, c.σ_major, c.coords)
end

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

fwhm_max(c::EllipticGaussianCovmat) = √(_eigen(c.covmat).values[2]) |> σ_to_fwhm
fwhm_min(c::EllipticGaussianCovmat) = √(_eigen(c.covmat).values[1]) |> σ_to_fwhm
fwhm_average(c::EllipticGaussianCovmat) = (_eigen(c.covmat).values |> prod)^0.25 |> σ_to_fwhm

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
    E = _eigen(c.covmat)
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
@accessor components(m::ModelComponent) = tuple(m)

Base.:(==)(a::MultiComponentModel, b::MultiComponentModel) = components(a) == components(b)

Base.:+(args::Union{ModelComponent,MultiComponentModel}...) = MultiComponentModel(reduce((a, b) -> (a..., b...), components.(args)))
Base.:*(x::Number, m::ModelComponent) = @modify(c -> x*c, flux(m))
Base.:*(x::Number, m::MultiComponentModel) = @modify(c -> x*c, components(m)[∗])

flux(m::MultiComponentModel) = sum(flux, components(m))

intensity_peak(c::Point) = flux(c) / effective_area(c)
intensity_peak(c::CircularGaussian) = flux(c) / effective_area(c)
intensity_peak(c::EllipticGaussian) = flux(c) / effective_area(c)
intensity_peak(c::EllipticGaussianCovmat) = flux(c) / effective_area(c)

intensity(c, uv::XYType) = intensity(c)(uv)

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

# XXX: the sum(...) method should be the only one, but it has performance issues for Tuples for some reason
function intensity(m::MultiComponentModel{<:Tuple})
    bycomp = map(intensity, components(m))
    (xy::XYType) -> map(f -> f(xy), bycomp) |> sum
end
function intensity(m::MultiComponentModel{<:AbstractArray})
    bycomp = map(intensity, components(m))
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

visibility(c, uv::UVType) = visibility(c)(uv)
visibility(c::Point) = @inline (uv::UVType) -> flux(c) * cis(visibility(angle, c, uv))
function visibility(c::CircularGaussian)
    mul = -2π^2 * c.σ^2
    @inline (uv::UVType) -> flux(c) * exp(mul * dot(uv, uv)) * cis(visibility(angle, c, uv))
end
visibility(c::EllipticGaussian) = visibility(EllipticGaussianCovmat(c))
visibility(c::EllipticGaussianCovmat) = @inline uv::UVType -> flux(c) * exp(-2π^2 * dot(uv, c.covmat, uv)) * cis(visibility(angle, c, uv))
visibility(m::MultiComponentModel{<:Tuple{Any}}) = visibility(only(components(m)))

# XXX: the sum(...) method should be the only one, but it has performance issues for Tuples for some reason
function visibility(m::MultiComponentModel{<:Tuple})
    bycomp = map(visibility, components(m))
    (xy::UVType) -> map(f -> f(xy), bycomp) |> sum
end
function visibility(m::MultiComponentModel{<:AbstractArray})
    bycomp = map(visibility, components(m))
    (xy::UVType) -> sum(f -> f(xy), bycomp)
end

visibility_envelope(::typeof(abs), c::Point, uvdist::Real) = c.flux ± eps(float(c.flux))
visibility_envelope(::typeof(abs), c::CircularGaussian, uvdist::Real) = c.flux * exp(-2π^2 * c.σ^2 * uvdist^2) ± eps(float(c.flux))
visibility_envelope(::typeof(abs), c::EllipticGaussian, uvdist::Real) = (c.flux * exp(-2π^2 * c.σ_major^2 * uvdist^2)) .. (c.flux * exp(-2π^2 * (c.σ_major * c.ratio_minor_major)^2 * uvdist^2))

visibility_envelope(f::typeof(abs), m::MultiComponentModel{<:Tuple{Any}}, uvdist::Real) = visibility_envelope(f, only(components(m)), uvdist)
visibility_envelope(f::typeof(angle), m::MultiComponentModel{<:Tuple{Any}}, uvdist::Real) = visibility_envelope(f, only(components(m)), uvdist)

function visibility_envelope(f::Function, model, uvdist::Real; npoints=100, seed=npoints)
	rng = copy(Random.default_rng())
	Random.seed!(rng, seed)
	ends = extrema(rand(rng, 0..2π, npoints)) do θ
		visibility(f, model, uvdist * SVector(sincos(θ)))
	end
    m = maximum(abs, ends)
    Interval(ends[1] - 0.01m, ends[2] + 0.01m)  # XXX: arbitrary margin, not guaranteed
end
function visibility_envelope(f::typeof(abs), model, uvdist::Real; kwargs...)
    i = @invoke visibility_envelope(f::Function, model, uvdist; kwargs...)
    @modify(e -> max(e, zero(e)), leftendpoint(i))
end
function visibility_envelope(f::typeof(angle), model, uvdist::Real; kwargs...)
    i = @invoke visibility_envelope(f::Function, model, uvdist; kwargs...)
    maxval = maximum(abs, endpoints(i))
    return -maxval .. maxval
end

visibility(f::Function, args...; kwargs...) = f(visibility(args...; kwargs...))
Broadcast.broadcasted(::typeof(visibility), f::Function, args...; kwargs...) = f.(visibility.(args...; kwargs...))


Unitful.ustrip(x::ModelComponent) = @modify(x -> ustrip.(x), x |> Properties())
Unitful.ustrip(x::MultiComponentModel) = @modify(ustrip, components(x)[∗])


Base.isapprox(a::ModelComponent, b::ModelComponent; kwargs...) =
    propertynames(a) == propertynames(b) && all(k -> isapprox(getproperty(a, k), getproperty(b, k); kwargs...), propertynames(a))
Base.isapprox(a::MultiComponentModel, b::MultiComponentModel; kwargs...) =
    length(components(a)) == length(components(b)) && all(((x, y),) -> isapprox(x, y; kwargs...), zip(components(a), components(b)))


AccessorsExtra.init_for_construct(::Type{Point}) = Point(NaN, SVector(NaN, NaN))
AccessorsExtra.init_for_construct(::Type{CircularGaussian}) = CircularGaussian(NaN, NaN, SVector(NaN, NaN))
AccessorsExtra.init_for_construct(::Type{EllipticGaussian}) = EllipticGaussian(NaN, NaN, NaN, NaN, SVector(NaN, NaN))

AccessorsExtra.@define_construct_by_set Point (typeof(flux), typeof(coords))
AccessorsExtra.@define_construct_by_set CircularGaussian (
    typeof(flux),
    Union{PropertyLens{:σ},typeof(fwhm_average),typeof(fwhm_max),typeof(fwhm_min),typeof(effective_area)},
    typeof(coords)
)
AccessorsExtra.@define_construct_by_set EllipticGaussian (
    typeof(flux),
    typeof(fwhm_max),
    Union{PropertyLens{:ratio_minor_major},typeof(fwhm_min)},
    typeof(position_angle),
    typeof(coords)
)

construct(::Type{CircularGaussian},
    p1::Pair{typeof(flux)},
    (o2, v2)::Pair{<:Union{typeof(intensity_peak),Base.Fix2{typeof(Tb_peak)}}},
    p3::Pair{typeof(coords)}
) = AccessorsExtra.construct_by_set(CircularGaussian, (p1, modifying(@o _.σ)(o2) => v2, p3))
