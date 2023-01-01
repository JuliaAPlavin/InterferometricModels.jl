using Test
using InterferometricModels
using Unitful, UnitfulAstro
using StaticArrays
using LinearAlgebra: normalize
using RectiGrids
import InterferometricModels as IM
import VLBIData as VLBI


@testset "components" begin
    @testset "circular" begin
        c = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
        @test intensity_peak(c) ≈ 1.5 / (2π*0.1^2)
        @test intensity(c)(SVector(1, 2)) ≈ intensity_peak(c)
        @test fwhm_max(c) ≈ fwhm_min(c) ≈ fwhm_average(c)
        w = fwhm_max(c)/2
        @test intensity(c)(SVector(1 + w, 2)) ≈ 0.5 * intensity_peak(c)
        @test intensity(c)(SVector(1, 2 - w)) ≈ 0.5 * intensity_peak(c)
        @test intensity(c)(SVector(1 + w/√(2), 2 + w/√(2))) ≈ 0.5 * intensity_peak(c)

        xs = range(-10, 10, length=2000)
        img = intensity(c).(grid(SVector, xs, xs))
        @test sum(img) * (20 / 2000)^2 ≈ 1.5  rtol=1e-3
        @test 0.99 <= maximum(img) / intensity_peak(c) <= 1
        @test SVector(map((a, i) -> a[i], axiskeys(img), Tuple(argmax(img)))) ≈ SVector(1, 2)  atol=0.03
    end

    @testset "elliptical" begin
        c = EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.))
        @test intensity_peak(c) ≈ 1.5 / (2π*0.5*0.25)
        @test intensity(c)(SVector(1, 2)) ≈ intensity_peak(c)
        @test fwhm_max(c) ≈ 0.5 * √(8 * log(2))
        @test fwhm_min(c) ≈ 0.25 * √(8 * log(2))
        @test fwhm_average(c) ≈ sqrt(0.5*0.25) * √(8 * log(2))
        
        off = normalize(SVector(0.3, 1)) * fwhm_max(c)/2
        @test intensity(c)(SVector(1, 2) + off) ≈ 0.5 * intensity_peak(c)
        off = normalize(SVector(1, -0.3)) * fwhm_min(c)/2
        @test intensity(c)(SVector(1, 2) + off) ≈ 0.5 * intensity_peak(c)

        xs = range(-10, 10, length=2000)
        img = intensity(c).(grid(SVector, xs, xs))
        @test sum(img) * (20 / 2000)^2 ≈ 1.5  rtol=1e-3
        @test 0.99 <= maximum(img) / intensity_peak(c) <= 1
        @test SVector(map((a, i) -> a[i], axiskeys(img), Tuple(argmax(img)))) ≈ SVector(1, 2)  atol=0.03
    end

    @testset "elliptical covmat" begin
        cs = [
            CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.)),
            EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.))
        ]
        for c in cs
            ccov = EllipticGaussianCovmat(c)
            @test flux(ccov) == flux(c)
            @test coords(ccov) == coords(c)
            @test intensity_peak(ccov) ≈ intensity_peak(c)
            xs = [[SVector(1., 2.), SVector(0., 0.)]; SVector.(randn(10), randn(10))]
            @test intensity(ccov).(xs) ≈ intensity(c).(xs)
        end
    end

    @testset "multicomponent" begin
        c1 = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
        c2 = EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.))
        cs = (c1, c1, c2)
        m = MultiComponentModel(cs)
        xs = [[SVector(1., 2.), SVector(0., 0.)]; SVector.(randn(10), randn(10))]
        @test intensity(m).(xs) ≈ 2 .* intensity(c1).(xs) .+ intensity(c2).(xs)
    end

    @testset "convolve beam" begin
        c = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
        @assert !InterferometricModels.is_beam(c)
        @assert InterferometricModels.is_beam(beam(CircularGaussian, σ=0.5))
        cc = convolve(c, beam(CircularGaussian, σ=0.5))
        @test cc isa CircularGaussian
        @test coords(cc) ≈ SVector(1, 2)
        @test 1.4 < intensity_peak(cc) < 1.5

        cc = convolve(c, beam(EllipticGaussian, σ_major=1., ratio_minor_major=0.5, pa_major=0.))
        @test cc isa EllipticGaussianCovmat
        @test coords(cc) ≈ SVector(1, 2)
        @test 1.4 < intensity_peak(cc) < 1.5

        c = EllipticGaussian(flux=1.5, σ_major=0.1, ratio_minor_major=0.7, pa_major=0.123, coords=SVector(1., 2.))
        cc = convolve(c, beam(CircularGaussian, σ=0.5))
        @test cc isa EllipticGaussianCovmat
        @test coords(cc) ≈ SVector(1, 2)
        @test 1.4 < intensity_peak(cc) < 1.5

        cc = convolve(c, beam(EllipticGaussian, σ_major=1., ratio_minor_major=0.5, pa_major=0.))
        @test cc isa EllipticGaussianCovmat
        @test coords(cc) ≈ SVector(1, 2)
        @test 1.4 < intensity_peak(cc) < 1.5
    end
end

# @testset begin
#     dmod = VLBI.load(VLBI.DifmapModel, "./data/difmap_model.mod")
#     mod = IM.model_from_difmap(dmod)
#     @test collect(map(flux, components(mod))) ≈ [0.649486, -0.00937331, 1.32054e9]
#     @test collect(map(intensity_peak, components(mod))) ≈ [Inf, -0.00022640511888964438, 0.011191592828378834]
#     @test collect(map(c -> Tb_peak(c, 15u"GHz"), components(mod))) ≈ [Inf, -1.3934202495658437e6, 6.887914967841405e7]

#     @test position_angle(components(mod)[1]) ≈ -0.8165976691060989
#     @test_throws MethodError position_angle(components(mod)[2])
#     @test position_angle(components(mod)[1] => components(mod)[2]) ≈ 0.510276455129927

#     @test_broken visibility(mod, SVector(0, 0))
#     @test_broken visibility(mod.components[1], SVector(0, 0))
#     @test visibility(mod.components[2], SVector(0, 0)) ≈ -0.00937331 + 0.0im
# end


import CompatHelperLocal as CHL
CHL.@check()
import Aqua
@testset "Aqua" begin
    Aqua.test_ambiguities(IM, recursive=false)
    Aqua.test_unbound_args(IM)
    Aqua.test_undefined_exports(IM)
    Aqua.test_stale_deps(IM)
end
