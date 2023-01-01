module InterferometricModels

using Unitful, UnitfulAstro
using IntervalSets
using StaticArrays
using LinearAlgebra
using Accessors

export
    ModelComponent, Point, CircularGaussian, EllipticGaussian, EllipticGaussianCovmat,
    MultiComponentModel, components,
    flux, coords, Tb_peak, fwhm_max, fwhm_average, fwhm_min, effective_area, position_angle, separation,
    intensity_peak, intensity,
    visibility, visibility_envelope,
    convolve, beam,
    ustrip

include("utils.jl")
include("models.jl")
include("convolve.jl")
include("from_vec.jl")

end
