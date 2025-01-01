struct Beam{C<:ModelComponent}
    comp::C
end

beam(T::Type; kwargs...) = Beam(T; kwargs...)

ModelComponent(b::Beam) = b.comp

Beam(::Type{CircularGaussian}; σ) = Beam(CircularGaussian(; flux=2π*σ^2, σ, coords=SVector(zero(σ), zero(σ))))
Beam(::Type{EllipticGaussian}; σ_major, ratio_minor_major, pa_major) = Beam(
    EllipticGaussian(; flux=2π*σ_major^2*ratio_minor_major, σ_major, ratio_minor_major, pa_major, coords=SVector(zero(σ_major), zero(σ_major)))
)

Base.getproperty(b::Beam, s::Symbol) = hasfield(typeof(b), s) ? getfield(b, s) : getproperty(b.comp, s)

Accessors.setproperties(b::Beam{<:CircularGaussian}, patch::NamedTuple{(:σ,)}) = Beam(normalize_as_beam(setproperties(b.comp, patch)))
Accessors.setproperties(b::Beam{<:EllipticGaussian}, patch::NamedTuple{(:ratio_minor_major,)}) = Beam(normalize_as_beam(setproperties(b.comp, patch)))

for f in (:fwhm_average, :fwhm_max, :fwhm_min, :effective_area, :position_angle)
    @eval $f(b::Beam) = $f(b.comp)
    @eval Accessors.set(b::Beam, ::typeof($f), val) = Beam(normalize_as_beam(set(b.comp, $f, val)))
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

normalize_as_beam(c::ModelComponent) = @p let
    c
    @set __.coords = SVector(zero(eltype(__.coords)), zero(eltype(__.coords)))
    set(__, modifying(@o _.flux)(intensity_peak), 1)
end
