module InterferometricModels

using Unitful, UnitfulAstro
using IntervalSets
using StaticArrays
using LinearAlgebra
using DataPipes
using AccessorsExtra
import AccessorsExtra: set, construct

export
    ModelComponent, Point, CircularGaussian, EllipticGaussian, EllipticGaussianCovmat,
    MultiComponentModel, components,
    flux, coords, Tb_peak, fwhm_max, fwhm_average, fwhm_min, effective_area, position_angle, separation,
    intensity_peak, intensity,
    visibility, visibility_envelope,
    convolve, beam, Beam,
    ustrip

include("utils.jl")
include("models.jl")
include("beam.jl")
include("convolve.jl")
include("modify.jl")


function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f === visibility_envelope && argtypes[1] == typeof(abs) && argtypes[2] <: MultiComponentModel
                println(io, "\n\nvisibility_envelope(::MultiComponentModel) requires IntervalArithmetic.jl: install and load that package")
            end
        end
    end
end

function _eigen(M::AbstractMatrix)
    Mu = unit(eltype(M))
    E_nou = eigen(ustrip(M))
    (values=E_nou.values*Mu, vectors=E_nou.vectors)
end

end
