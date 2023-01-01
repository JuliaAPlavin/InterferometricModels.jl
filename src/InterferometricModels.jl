module InterferometricModels

using Unitful, UnitfulAstro
using IntervalSets
using StaticArrays
using LinearAlgebra: det, dot, Diagonal

export
    ModelComponent, Point, CircularGaussian, EllipticGaussian, EllipticGaussianCovmat,
    flux, coords, Tb_peak, fwhm_max, fwhm_average, fwhm_min, effective_area, intensity_peak, intensity, position_angle,
    visibility, visibility_amplitude, visibility_phase, visibility_envelope,
    MultiComponentModel, components

include("utils.jl")
include("models.jl")
include("from_vec.jl")

end
