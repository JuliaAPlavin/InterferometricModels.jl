using Test
using InterferometricModels
using Unitful, UnitfulAstro
using StaticArrays
using LinearAlgebra
using RectiGrids
using IntervalSets


@testset "point" begin
    c = Point(flux=1.5, coords=SVector(1., 2.))

    @test fwhm_max(c) == fwhm_min(c) == fwhm_average(c) == 0
    @test intensity_peak(c) ≈ Inf
    @test_broken intensity(c)(SVector(1, 2)) ≈ Inf
    @test_broken intensity(c)(SVector(1.1, 2)) == 0

    @test visibility(c, SVector(0, 0)) == flux(c)
    @test visibility(c, SVector(-1.23, 4.56)) ≈ flux(c) * cis(angle(0.01415 - 0.01170im))  rtol=1e-3
    @test visibility(abs, c, SVector(-1.23, 4.56)) == flux(c)
    @test mod2pi(visibility(angle, c, SVector(-1.23, 4.56))) ≈ mod2pi(angle(0.01415 - 0.01170im))  rtol=1e-4
end

@testset "circular" begin
    c = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))

    @test intensity_peak(c) ≈ 1.5 / (2π*0.1^2)
    @test intensity(c)(SVector(1, 2)) ≈ intensity_peak(c)
    @test fwhm_max(c) ≈ fwhm_min(c) ≈ fwhm_average(c)
    w = fwhm_max(c)/2
    @test intensity(c)(SVector(1 + w, 2)) ≈ 0.5 * intensity_peak(c)
    @test intensity(c)(SVector(1, 2 - w)) ≈ 0.5 * intensity_peak(c)
    @test intensity(c)(SVector(1 + w/√(2), 2 + w/√(2))) ≈ 0.5 * intensity_peak(c)

    @test visibility(c, SVector(0, 0)) == flux(c)
    @test visibility(c, SVector(1, 1/2)) ≈ 1.1720  rtol=1e-4
    @test visibility(c, SVector(-1.23, 4.56)) ≈ 0.0141454 - 0.0117022im  rtol=1e-4
    @test visibility(abs, c, SVector(-1.23, 4.56)) ≈ abs(0.0141454 - 0.0117022im)  rtol=1e-4
    @test mod2pi(visibility(angle, c, SVector(-1.23, 4.56))) ≈ mod2pi(angle(0.01415 - 0.01170im))  rtol=1e-4

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

    @test visibility(c, SVector(0, 0)) == flux(c)
    @test visibility(c, SVector(1, 1/2)) ≈ 0.036524  rtol=1e-4
    @test visibility(c, SVector(-0.123, 0.456)) ≈ 0.15221 - 0.60867im  rtol=1e-4
    @test visibility(abs, c, SVector(-0.123, 0.456)) ≈ abs(0.15221 - 0.60867im)  rtol=1e-4
    @test mod2pi(visibility(angle, c, SVector(-0.123, 0.456))) ≈ mod2pi(angle(0.15221 - 0.60867im))  rtol=1e-4

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
        @test visibility.(ccov, xs) ≈ visibility.(c, xs)
        @test visibility.(abs, ccov, xs) ≈ visibility.(abs, c, xs)
        @test visibility.(angle, ccov, xs) ≈ visibility.(angle, c, xs)
    end
end

@testset "multicomponent" begin
    c1 = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
    c2 = EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(0.5, -0.1))
    @test position_angle(c1 => c2) ≈ -2.907849472720892
    @test position_angle(c2 => c1) ≈ 0.23374318086890142
    @test separation(c1, c2) ≈ 2.1587033144922905
    cs = (c1, c1, c2)
    m = MultiComponentModel(cs)
    xs = [[SVector(1., 2.), SVector(0., 0.)]; SVector.(randn(10), randn(10))]
    @test intensity(m).(xs) ≈ 2 .* intensity(c1).(xs) .+ intensity(c2).(xs)
    @test visibility.(m, xs) ≈ 2 .* visibility.(c1, xs) .+ visibility.(c2, xs)
    @test flux(m) == 4.5
end

@testset "convolve beam" begin
    c = Point(flux=1.5, coords=SVector(1., 2.))
    cc = convolve(c, beam(CircularGaussian, σ=0.5))
    @test cc isa CircularGaussian
    @test coords(cc) ≈ SVector(1, 2)
    @test intensity_peak(cc) == 1.5

    cc = convolve(c, beam(EllipticGaussian, σ_major=1., ratio_minor_major=0.5, pa_major=0.))
    @test cc isa EllipticGaussian
    @test coords(cc) ≈ SVector(1, 2)
    @test intensity_peak(cc) == 1.5

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

    cm = convolve(MultiComponentModel((c,)), beam(CircularGaussian, σ=0.5))
    @test only(components(cm)) == convolve(c, beam(CircularGaussian, σ=0.5))
end

@testset "envelopes" begin
    cs = [
        Point(flux=1.5, coords=SVector(1., 2.)),
        CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.)),
        EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.)),
    ]
    # append!(cs, EllipticGaussianCovmat.(cs[2:end]))
    xs = [[SVector(1., 2.), SVector(0., 0.)]; SVector.(randn(10), randn(10))]
    @testset for c in cs, x in xs
        @test visibility(abs, c, x) ∈ visibility_envelope(abs, c, norm(x))
        @test mod2pi(visibility(angle, c, x)+π)-π ∈ visibility_envelope(angle, c, norm(x))
    end
end

@testset "ustrip" begin
    c = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
    cu = CircularGaussian(flux=1.5u"Jy", σ=0.1u"°", coords=SVector(1., 2.)*u"°")
    @test ustrip(c) == c
    @test ustrip(cu) == c
    m = MultiComponentModel((c,))
    mu = MultiComponentModel((cu,))
    @test ustrip(m) == m
    @test ustrip(mu) == m

    @test ustrip(1u"m"..2u"m") == 1..2
    @test ustrip(u"cm", 1u"m"..2u"m") == 100..200
end


import CompatHelperLocal as CHL
CHL.@check()
import Aqua
@testset "Aqua" begin
    Aqua.test_ambiguities(InterferometricModels, recursive=false)
    Aqua.test_unbound_args(InterferometricModels)
    Aqua.test_undefined_exports(InterferometricModels)
    Aqua.test_stale_deps(InterferometricModels)
end
