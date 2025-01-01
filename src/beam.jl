struct Beam{C<:ModelComponent}
    comp::C
end

beam(T::Type; kwargs...) = Beam(T; kwargs...)

ModelComponent(b::Beam) = b.comp

Beam(::Type{CircularGaussian}; σ, keep_flux=false) = Beam(normalize_as_beam(CircularGaussian(; flux=true, σ, coords=SVector(zero(σ), zero(σ))); keep_flux))
Beam(::Type{EllipticGaussian}; σ_major, ratio_minor_major, pa_major, keep_flux=false) = Beam(normalize_as_beam(
    EllipticGaussian(; flux=true, σ_major, ratio_minor_major, pa_major, coords=SVector(zero(σ_major), zero(σ_major)));
    keep_flux))

Base.getproperty(b::Beam, s::Symbol) = hasfield(typeof(b), s) ? getfield(b, s) : getproperty(b.comp, s)

Accessors.setproperties(b::Beam{<:CircularGaussian}, patch::NamedTuple{(:σ,)}) = Beam(normalize_as_beam(setproperties(b.comp, patch); keep_flux=_keep_flux(b)))
Accessors.setproperties(b::Beam{<:EllipticGaussian}, patch::NamedTuple{(:ratio_minor_major,)}) = Beam(normalize_as_beam(setproperties(b.comp, patch); keep_flux=_keep_flux(b)))

for f in (:flux, :fwhm_average, :fwhm_max, :fwhm_min, :effective_area, :position_angle)
    @eval $f(b::Beam) = $f(b.comp)
    @eval Accessors.set(b::Beam, ::typeof($f), val) = Beam(normalize_as_beam(set(b.comp, $f, val); keep_flux=_keep_flux(b)))
end

AccessorsExtra.init_for_construct(::Type{Beam{C}}) where {C<:ModelComponent} = Beam(AccessorsExtra.init_for_construct(C))

AccessorsExtra.@define_construct_by_set Beam{CircularGaussian} (
    Union{PropertyLens{:σ},typeof(fwhm_average),typeof(fwhm_max),typeof(fwhm_min),typeof(effective_area)},
)
AccessorsExtra.@define_construct_by_set Beam{EllipticGaussian} (
    typeof(fwhm_max),
    Union{PropertyLens{:ratio_minor_major},typeof(fwhm_min)},
    typeof(position_angle),
)

convolve(model, b::Beam) = convolve(model, b.comp)

is_beam(c::ModelComponent) = false
is_beam(b::Beam) = true
is_flux_unity(b::Beam) = is_flux_unity(b.comp)
is_flux_unity(c::ModelComponent) = flux(c) ≈ one(flux(c))
is_intensity_unity(b::Beam) = is_intensity_unity(b.comp)
is_intensity_unity(c::ModelComponent) = intensity_peak(c) ≈ one(intensity_peak(c))
function _keep_flux(b::Beam)
    is_flux_unity(b) && return true
    is_intensity_unity(b) && return false
    isnan(flux(b)) && return false
    error("Beam should have either flux of unity (got $(flux(b))) or intensity peak of unity (got $(intensity_peak(b.comp)))")
end

normalize_as_beam(c::ModelComponent; keep_flux=false) = @p let
    c
    @set __.coords = SVector(zero(eltype(__.coords)), zero(eltype(__.coords)))
    keep_flux ? (@set __.flux = true) :
        set(__, modifying(@o _.flux)(intensity_peak), 1)
end

Unitful.ustrip(x::Beam) = @modify(ustrip, x.comp)
