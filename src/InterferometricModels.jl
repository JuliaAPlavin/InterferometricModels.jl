module InterferometricModels

using Unitful, UnitfulAstro
using IntervalSets
using StaticArrays
using LinearAlgebra
using DataPipes
using AccessorsExtra
import AccessorsExtra: set, construct
using AccessorsExtra.InverseFunctions: square
using StructHelpers
import Random
using AccessorsExtra: ConstrainedLens, FixArgs

export
    ModelComponent, Point, CircularGaussian, EllipticGaussian, EllipticGaussianCovmat,
    MultiComponentModel, components,
    flux, coords, Tb_peak, fwhm_max, fwhm_average, fwhm_min, effective_area, position_angle, separation,
    intensity_peak, intensity,
    visibility, visibility_envelope,
    convolve, beam, Beam,
    ustrip,
    UV, UVW

include("utils.jl")
include("models.jl")
include("beam.jl")
include("convolve.jl")
include("modify.jl")

function _eigen(M::AbstractMatrix)
    Mu = unit(eltype(M))
    E_nou = eigen(ustrip(M))
    (values=E_nou.values*Mu, vectors=E_nou.vectors)
end

end
