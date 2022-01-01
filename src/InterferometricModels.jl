module InterferometricModels

using Parameters
using Unitful, UnitfulAstro
using IntervalSets
using StaticArrays
using Utils  # for u"mas"

export
    ModelComponent, PointSource, CircularGaussian, EllipticGaussian,
    flux, coords, Tb_peak, fwhm_max, fwhm_average, fwhm_min, effective_area, intensity_peak, intensity,
    visibility, visibility_amplitude, visibility_phase, visibility_envelope,
    MultiComponentModel, components

include("utils.jl")
include("models.jl")

end
