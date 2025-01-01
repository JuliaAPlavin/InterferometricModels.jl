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
include("modify.jl")


function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            # @info "" exc.f argtypes kwargs
            if exc.f === visibility_envelope && argtypes[1] == typeof(abs) && argtypes[2] <: MultiComponentModel
                println(io, "\n\nvisibility_envelope(::MultiComponentModel) requires IntervalArithmetic.jl: install and load that package")
            end
        end
    end
end

end
