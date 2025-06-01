using TestItems
using TestItemRunner
@run_package_tests

@testitem "point" begin
    using StaticArrays
    using AccessorsExtra: construct, test_construct_laws
    using Unitful

    c = Point(flux=1.5, coords=SVector(1., 2.))
    @test c === Point(flux=1.5, coords=(1., 2.))
    @test c === Point(flux=1.5, coords=[1., 2.])
    @test c === typeof(c)(flux=1.5, coords=SVector(1., 2.))

    @test c ≈ Point(flux=1.5f0 * 1.00001f0, coords=SVector(1f0, 2f0) * 1.00001f0)
    @test !(c ≈ Point(flux=1.6, coords=SVector(1f0, 2f0)))

    @test fwhm_max(c) == fwhm_min(c) == fwhm_average(c) == 0
    @test intensity_peak(c) ≈ Inf
    @test_broken intensity(c)(SVector(1, 2)) ≈ Inf
    @test_broken intensity(c)(SVector(1.1, 2)) == 0

    @test visibility(c, SVector(0, 0)) == flux(c)
    @test visibility(c, SVector(-1.23, 4.56)) ≈ flux(c) * cis(angle(0.01415 - 0.01170im))  rtol=1e-3
    @test visibility(abs, c, SVector(-1.23, 4.56)) == flux(c)
    
    @test visibility.(c, [SVector(0, 0)]) |> only == flux(c)
    @test visibility.(c, [SVector(-1.23, 4.56)]) |> only ≈ flux(c) * cis(angle(0.01415 - 0.01170im))  rtol=1e-3
    @test visibility.(abs, c, [SVector(-1.23, 4.56)]) |> only == flux(c)

    @test mod2pi(visibility(angle, c, SVector(-1.23, 4.56))) ≈ mod2pi(angle(0.01415 - 0.01170im))  rtol=1e-4

    test_construct_laws(Point, flux=>2.5, coords=>SVector(1, 2))
    test_construct_laws(Point, flux=>2.5u"m", coords=>SVector(1, 2))
    test_construct_laws(Point, flux=>2.5, coords=>SVector(1, 2)u"°")
    test_construct_laws(Point, flux=>2.5u"m", coords=>SVector(1, 2)u"°")
end

@testitem "circular" begin
    using StaticArrays
    using RectiGrids
    using AccessorsExtra: @o, construct, test_construct_laws
    using Unitful

    c = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
    @test c === CircularGaussian(flux=1.5, σ=0.1, coords=(1., 2.))
    @test c == CircularGaussian(flux=1.5, σ=0.1, coords=[1, 2])

    @test c ≈ CircularGaussian(flux=1.5f0 * 1.00001f0, σ=0.1f0 * 1.00001f0, coords=SVector(1f0, 2f0))
    @test !(c ≈ CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1.3, 2.)))

    @test CircularGaussian(Point(flux=1.5f0, coords=SVector(1., 2.))) === CircularGaussian(flux=1.5f0, σ=0., coords=SVector(1., 2.))
    @test CircularGaussian(Point(flux=1.5f0, coords=SVector(1f0, 2f0))) === CircularGaussian(flux=1.5f0, σ=0f0, coords=SVector(1f0, 2f0))
    for c in Any[Point(flux=1.5, coords=SVector(1., 2.)),
                 Point(flux=1.5f0, coords=SVector(1., 2.)),
                 Point(flux=1.5f0, coords=SVector(1f0, 2f0))]
        @test (@inferred Point(CircularGaussian(c))) === c
    end

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

    @test construct(CircularGaussian, flux=>2.5, fwhm_average=>0.2, coords=>SVector(1, 2)) ≈
          CircularGaussian(flux=2.5, σ=InterferometricModels.fwhm_to_σ(0.2), coords=SVector(1, 2))
    test_construct_laws(CircularGaussian, flux=>2.5, fwhm_average=>0.2, coords=>SVector(1, 2))
    test_construct_laws(CircularGaussian, flux=>2.5, fwhm_max=>0.2u"°", coords=>SVector(1, 2)u"°")
    test_construct_laws(CircularGaussian, flux=>2.5u"m", fwhm_average=>0.2, coords=>SVector(1, 2))
    test_construct_laws(CircularGaussian, flux=>2.5u"m", (@o _.σ)=>0.2u"°", coords=>SVector(1, 2))
    test_construct_laws(CircularGaussian, flux=>2.5u"m", effective_area=>0.2u"°^2", coords=>SVector(1, 2); cmp=(≈))
    test_construct_laws(CircularGaussian, flux=>2.5u"m", intensity_peak=>0.2u"m/°^2", coords=>SVector(1, 2); cmp=(≈))
    test_construct_laws(CircularGaussian, flux=>2.5u"W/m^2/Hz", (@o Tb_peak(_, 5u"GHz"))=>0.2u"K", coords=>SVector(1, 2); cmp=(≈))
end

@testitem "elliptical" begin
    using StaticArrays
    using RectiGrids
    using LinearAlgebra: normalize
    using AccessorsExtra: @o, construct, test_construct_laws
    using Unitful

    c = EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.))
    @test c === EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=(1., 2.))
    @test c == EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=[1, 2])
    @test c == EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=1//2, pa_major=deg2rad(16.6992), coords=(1, 2))

    @test c ≈ EllipticGaussian(flux=1.5f0 * 1.00001f0, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.))
    @test c ≈ EllipticGaussian(flux=1.5f0, σ_major=0.5f0, ratio_minor_major=0.5f0, pa_major=Float32(deg2rad(16.6992f0)), coords=SVector(1f0, 2f0))
    @test c ≈ EllipticGaussian(flux=1.5f0, σ_major=0.25f0, ratio_minor_major=2f0, pa_major=deg2rad(-90+16.6992f0), coords=SVector(1f0, 2f0))
    @test c ≈ EllipticGaussian(flux=1.5, σ_major=0.25f0, ratio_minor_major=2f0, pa_major=deg2rad(-90+16.6992f0), coords=SVector(1, 2))
    @test !(c ≈ EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.9), coords=SVector(1., 2.)))

    @test (@inferred EllipticGaussian(Point(flux=1.5f0, coords=SVector(1., 2.)))) === EllipticGaussian(flux=1.5f0, σ_major=0., ratio_minor_major=1., pa_major=0., coords=SVector(1., 2.))
    @test (@inferred EllipticGaussian(Point(flux=1.5f0, coords=SVector(1f0, 2f0)))) === EllipticGaussian(flux=1.5f0, σ_major=0f0, ratio_minor_major=1f0, pa_major=0f0, coords=SVector(1f0, 2f0))
    @test (@inferred EllipticGaussian(CircularGaussian(flux=1.5f0, σ=0.1, coords=SVector(1., 2.)))) === EllipticGaussian(flux=1.5f0, σ_major=0.1, ratio_minor_major=1., pa_major=0., coords=SVector(1., 2.))
    @test (@inferred EllipticGaussian(CircularGaussian(flux=1.5f0, σ=0.1f0, coords=SVector(1f0, 2f0)))) === EllipticGaussian(flux=1.5f0, σ_major=0.1f0, ratio_minor_major=1f0, pa_major=0f0, coords=SVector(1f0, 2f0))
    for c in Any[Point(flux=1.5, coords=SVector(1., 2.)),
                 Point(flux=1.5f0, coords=SVector(1., 2.)),
                 Point(flux=1.5f0, coords=SVector(1f0, 2f0)),]
        @test (@inferred Point(EllipticGaussian(c))) === c
    end
    for c in Any[CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.)),
                 CircularGaussian(flux=1.5f0, σ=0.1f0, coords=SVector(1f0, 2f0)),]
        @test (@inferred CircularGaussian(EllipticGaussian(c))) === c
    end

    @test intensity_peak(c) ≈ 1.5 / (2π*0.5*0.25)
    @test intensity(c)(SVector(1, 2)) ≈ intensity_peak(c)
    @test fwhm_max(c) ≈ 0.5 * √(8 * log(2))
    @test fwhm_min(c) ≈ 0.25 * √(8 * log(2))
    @test fwhm_average(c) ≈ sqrt(0.5*0.25) * √(8 * log(2))

    c2 = EllipticGaussian(flux=1.5, σ_major=0.25, ratio_minor_major=2.0, pa_major=deg2rad(16.6992 + 90), coords=SVector(1., 2.))
    @test_broken EllipticGaussian(flux=1.5, σ_major=1, ratio_minor_major=2.0, pa_major=deg2rad(16.6992 + 90), coords=SVector(1., 2.))
    @test intensity_peak(c2) ≈ intensity_peak(c)
    @test fwhm_min(c2) ≈ fwhm_min(c)
    @test fwhm_max(c2) ≈ fwhm_max(c)
    @test fwhm_average(c2) ≈ fwhm_average(c)
    @test mod(position_angle(c2), π) ≈ mod(position_angle(c), π)
    @test EllipticGaussianCovmat(c2) ≈ EllipticGaussianCovmat(c)
    
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

    @test construct(EllipticGaussian, flux=>2.5, fwhm_max=>0.2, fwhm_min=>0.1, position_angle=>-0.5, coords=>SVector(1, 2)) ≈
          EllipticGaussian(flux=2.5, σ_major=InterferometricModels.fwhm_to_σ(0.2), ratio_minor_major=0.5, pa_major=-0.5, coords=SVector(1, 2))
    test_construct_laws(EllipticGaussian, flux=>2.5, fwhm_max=>0.2, (@o _.ratio_minor_major) => 0.5, position_angle=>-0.5, coords=>SVector(1, 2))
    test_construct_laws(EllipticGaussian, flux=>2.5, fwhm_max=>0.2, fwhm_min=>0.1, position_angle=>-0.5, coords=>SVector(1, 2))
    test_construct_laws(EllipticGaussian, flux=>2.5, fwhm_max=>0.2u"°", fwhm_min=>0.1u"°", position_angle=>-0.5, coords=>SVector(1, 2)u"°")
    test_construct_laws(EllipticGaussian, flux=>2.5u"m", fwhm_max=>0.2, fwhm_min=>0.1, position_angle=>-0.5, coords=>SVector(1, 2))
    test_construct_laws(EllipticGaussian, flux=>2.5u"m", fwhm_max=>0.2u"°", fwhm_min=>0.1u"°", position_angle=>-0.5, coords=>SVector(1, 2))
end

@testitem "elliptical covmat" begin
    using StaticArrays
    using Unitful

    cs = [
        CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.)),
        EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.)),
        EllipticGaussian(flux=1.5, σ_major=0.5u"rad", ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.)u"rad"),
        EllipticGaussian(flux=1.5, σ_major=0.5u"rad", ratio_minor_major=0.5, pa_major=16.6992u"°", coords=SVector(1., 2.)u"rad")
    ]
    @testset for c in cs
        ccov = EllipticGaussianCovmat(c)
        ccov_el = EllipticGaussian(ccov)
        @test flux(ccov_el) == flux(ccov) == flux(c)
        @test coords(ccov_el) == coords(ccov) == coords(c)
        @test fwhm_average(ccov_el) ≈ fwhm_average(ccov) ≈ fwhm_average(c)
        @test fwhm_max(ccov_el) ≈ fwhm_max(ccov) ≈ fwhm_max(c)
        @test fwhm_min(ccov_el) ≈ fwhm_min(ccov) ≈ fwhm_min(c)
        @test intensity_peak(ccov_el) ≈ intensity_peak(ccov) ≈ intensity_peak(c)
        xs = [[SVector(1., 2.), SVector(0., 0.)]; SVector.(randn(10), randn(10))]
        @test intensity(ccov_el).(xs) ≈ intensity(ccov).(xs) ≈ intensity(c).(xs)
        @test visibility.(ccov_el, xs) ≈ visibility.(ccov, xs) ≈ visibility.(c, xs)
        @test visibility.(abs, ccov, xs) ≈ visibility.(abs, c, xs)
        @test visibility.(angle, ccov, xs) ≈ visibility.(angle, c, xs)

        el = EllipticGaussian(c)
        @test ccov_el.σ_major ≈ el.σ_major
        @test ccov_el.ratio_minor_major ≈ el.ratio_minor_major
        @test ccov_el.pa_major ≈ el.pa_major
    end
end

@testitem "multicomponent" begin
    using StaticArrays
    using Accessors

    c1 = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
    c2 = EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(0.5, -0.1))
    @test position_angle(c1 => c2) ≈ -2.907849472720892
    @test position_angle(c2 => c1) ≈ 0.23374318086890142
    @test separation(c1, c2) ≈ 2.1587033144922905
    cs = (c1, c1, c2)
    @testset for m in (MultiComponentModel(cs), MultiComponentModel(collect(cs)))
        @test deepcopy(m) == m
        xs = [[SVector(1., 2.), SVector(0., 0.)]; SVector.(randn(10), randn(10))]
        @test intensity(m).(xs) ≈ 2 .* intensity(c1).(xs) .+ intensity(c2).(xs)
        @test visibility.(m, xs) ≈ 2 .* visibility.(c1, xs) .+ visibility.(c2, xs)
        @test visibility.(m, UV.(xs)) ≈ 2 .* visibility.(c1, UV.(xs)) .+ visibility.(c2, xs)
        @test flux(m) == 4.5
        @test deepcopy(m) == m
        @test @set(components(m)[2].flux = 1.500002f0) ≈ m
        @test !(@set(components(m)[2].flux = 1.7f0) ≈ m)
    end
end

@testitem "model arithmetic" begin
    using StaticArrays

    c1 = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
    c2 = EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(0.5, -0.1))
    c3 = CircularGaussian(flux=1, σ=0.1, coords=SVector(1, 2))
    @test c1 + c2 === MultiComponentModel((c1, c2))
    @test c1 + c2 + c3 === MultiComponentModel((c1, c2, c3))
    @test (c1 + c2) + c3 === MultiComponentModel((c1, c2, c3))
    @test c1 + (c2 + c3) === MultiComponentModel((c1, c2, c3))

    @test 2 * c1 === CircularGaussian(flux=3., σ=0.1, coords=SVector(1., 2.))
    @test 3 * (c1 + (c2 + c3)) === MultiComponentModel((3c1, 3c2, 3c3))
end

@testitem "beam" begin
    using AccessorsExtra
    using AccessorsExtra: test_construct_laws
    using Unitful, UnitfulAngles

    b = Beam(CircularGaussian, σ=0.5)
    @test b == construct(Beam{CircularGaussian}, fwhm_average => InterferometricModels.σ_to_fwhm(0.5))

    b = Beam(EllipticGaussian, σ_major=0.5, ratio_minor_major=0.5, pa_major=0.123)
    @test b == construct(Beam{EllipticGaussian}, fwhm_max => InterferometricModels.σ_to_fwhm(0.5), (@o _.ratio_minor_major) => 0.5, position_angle => 0.123)
    
    @test ustrip(b) == Beam(EllipticGaussian, σ_major=0.5, ratio_minor_major=0.5, pa_major=0.123)
    @test ustrip(Beam(EllipticGaussian, σ_major=0.5u"mas", ratio_minor_major=0.5, pa_major=0.123u"rad")) == Beam(EllipticGaussian, σ_major=0.5, ratio_minor_major=0.5, pa_major=0.123)

    for f in ((@o _.σ), fwhm_average, fwhm_max, fwhm_min)
        test_construct_laws(Beam{CircularGaussian}, f => 0.5, type=Beam{<:CircularGaussian})
        test_construct_laws(Beam{CircularGaussian}, f => 0.5u"mas", type=Beam{<:CircularGaussian})
    end
    test_construct_laws(Beam{CircularGaussian}, effective_area => 0.5, cmp=(≈), type=Beam{<:CircularGaussian})
    test_construct_laws(Beam{CircularGaussian}, effective_area => 0.5u"mas"^2, cmp=(≈), type=Beam{<:CircularGaussian})

    test_construct_laws(Beam{EllipticGaussian}, fwhm_max => 0.5, (@o _.ratio_minor_major) => 0.6, position_angle => 0.123, type=Beam{<:EllipticGaussian})
    test_construct_laws(Beam{EllipticGaussian}, fwhm_max => 0.5u"mas", (@o _.ratio_minor_major) => 0.6, position_angle => 0.123, type=Beam{<:EllipticGaussian})
    test_construct_laws(Beam{EllipticGaussian}, fwhm_max => 0.5u"mas", fwhm_min => 0.3u"mas", position_angle => 0.123, type=Beam{<:EllipticGaussian})
end

@testitem "convolve beam" begin
    using StaticArrays

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

@testitem "envelopes" begin
    using StaticArrays
    using LinearAlgebra: norm
    using Unitful
    using IntervalSets

    cs = [
        Point(flux=1.5, coords=SVector(1., 2.)),
        CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.)),
        EllipticGaussian(flux=1.5, σ_major=0.5, ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.)),
        MultiComponentModel((
            CircularGaussian(1, 0.1, SVector(0, 0)),
        )),
        MultiComponentModel((
            CircularGaussian(1, 0.1, SVector(0, 0)),
            CircularGaussian(1, 0.15, SVector(0.7, 0)),
        )),

        EllipticGaussian(flux=1.5u"W", σ_major=0.5u"rad", ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.)u"rad"),
        MultiComponentModel((
            CircularGaussian(1u"W", 0.1u"rad", SVector(0, 0)u"rad"),
            CircularGaussian(1u"W", 0.15u"rad", SVector(0.7, 0)u"rad"),
        )),
        MultiComponentModel((
            CircularGaussian(1, 0.1u"°", SVector(0, 0)u"°"),
            CircularGaussian(1, 0.15u"°", SVector(0.7, 0)u"°"),
        )),
    ]
    append!(cs, EllipticGaussianCovmat.(cs[2:3]))

    xs = [[SVector(1., 2.), SVector(0., 0.)]; SVector.(randn(10), randn(10))]
    @testset for c in cs, x in xs
        ival = visibility_envelope(abs, c, norm(x))
        @test visibility(abs, c, x) ∈ ival
        @test minimum(ival) ≥ zero(minimum(ival))

        ival = visibility_envelope(angle, c, norm(x))
        @test rem(visibility(angle, c, x), 2π, RoundNearest) ∈ ival
        @test minimum(ival) ≈ -maximum(ival)

        @test rem(visibility(rad2deg∘angle, c, x), 360, RoundNearest) ∈ visibility_envelope(rad2deg∘angle, c, norm(x))
        @test mod2pi(visibility(u"rad"∘angle, c, x)+π)-π ∈ visibility_envelope(u"rad"∘angle, c, norm(x))
        @test mod(visibility(u"°"∘angle, c, x)+180, 0..360)-180 ∈ visibility_envelope(u"°"∘angle, c, norm(x))
    end
end

@testitem "utils" begin
    using Accessors
    using Accessors.InverseFunctions
    using Unitful
    using UnitfulAstro
    using UnitfulAngles

    @test InterferometricModels.σ_to_fwhm(1)::Float64 ≈ 2.3548200450309493
    @test InterferometricModels.σ_to_fwhm(0.1)::Float64 ≈ 0.23548200450309495
    @test InterferometricModels.σ_to_fwhm(0.1f0)::Float32 ≈ 0.23548200450309495
    @test ustrip(InterferometricModels.σ_to_fwhm(0.1f0u"m"))::Float32 ≈ 0.23548200450309495
    InverseFunctions.test_inverse(InterferometricModels.σ_to_fwhm, 0.1)
    InverseFunctions.test_inverse(InterferometricModels.σ_to_fwhm, 0.1u"m")
    InverseFunctions.test_inverse(InterferometricModels.fwhm_to_σ, 0.1)
    InverseFunctions.test_inverse(InterferometricModels.fwhm_to_σ, 0.1u"m")
    
    f = @optic(InterferometricModels.intensity_to_Tb(_, 5u"GHz"))
    @test f(0.1u"Jy/mas^2") ≈ 5.539089534545945e9u"K"
    @test f(0.1) ≈ 5.539089534545945e9
    InverseFunctions.test_inverse(f, 0.1u"Jy/mas^2")
    InverseFunctions.test_inverse(f, 0.1)
end

@testitem "ustrip" begin
    using StaticArrays
    using Unitful, UnitfulAstro
    using IntervalSets

    c = CircularGaussian(flux=1.5, σ=0.1, coords=SVector(1., 2.))
    cu = CircularGaussian(flux=1.5u"Jy", σ=0.1u"°", coords=SVector(1., 2.)*u"°")
    @test ustrip(c) == c
    @test ustrip(cu) == c
    m = MultiComponentModel((c,))
    mu = MultiComponentModel((cu,))
    @test ustrip(m) == m
    @test ustrip(mu) == m
end

@testitem "set" begin
    using StaticArrays
    using Unitful, UnitfulAstro, UnitfulAngles
    using Accessors

    c = Point(flux=1.5u"Jy", coords=SVector(1., 2.)u"mas")
    @test flux(@set(flux(c) = 2u"Jy")) == 2u"Jy"
    @test coords(@set(coords(c) = SVector(-10, 0.5)u"mas")) == SVector(-10, 0.5)u"mas"

    c = CircularGaussian(flux=1.0u"Jy", σ=0.1u"mas", coords=SVector(0, 0.)u"mas")
    @test flux(@set(flux(c) = 2u"Jy")) == 2u"Jy"
    @test fwhm_max(@set(fwhm_max(c) = 1.0u"mas")) == 1u"mas"
    @test fwhm_min(@set(fwhm_min(c) = 1.0u"mas")) == 1u"mas"
    @test fwhm_average(@set(fwhm_average(c) = 1.0u"mas")) == 1u"mas"
    @test coords(@set(coords(c) = SVector(-10, 0.5)u"mas")) == SVector(-10, 0.5)u"mas"

    c = EllipticGaussian(flux=1.5u"Jy", σ_major=0.5u"mas", ratio_minor_major=0.5, pa_major=deg2rad(16.6992), coords=SVector(1., 2.)u"mas")
    @test flux(@set(flux(c) = 2u"Jy")) == 2u"Jy"
    @test fwhm_max(@set(fwhm_max(c) = 1.0u"mas")) == 1u"mas"
    @test fwhm_min(@set(fwhm_max(c) = 1.0u"mas")) == fwhm_min(c)
    @test fwhm_min(@set(fwhm_min(c) = 1.0u"mas")) == 1u"mas"
    @test fwhm_max(@set(fwhm_min(c) = 1.0u"mas")) == fwhm_max(c)
    @test coords(@set(coords(c) = SVector(-10, 0.5)u"mas")) == SVector(-10, 0.5)u"mas"
end

@testitem "set so that" begin
    using StaticArrays
    using Unitful, UnitfulAstro, UnitfulAngles
    using AccessorsExtra

    @testset for c0 in [
        CircularGaussian(flux=1.0u"Jy", σ=0.1u"mas", coords=SVector(0, 0.)u"mas"),
        CircularGaussian(flux=1.0u"Jy", σ=0.0u"mas", coords=SVector(0, 0.)u"mas"),
        CircularGaussian(flux=1.0u"Jy", σ=NaN*u"mas", coords=SVector(0, 0.)u"mas"),
    ]
        @testset "$o $func" for
                (o, o_other) in [
                    (@optic(_.σ), @optic(_.flux)),
                    (@optic(_.flux), @optic(_.σ)),
                ],
                (func, tgt) in [
                    intensity_peak => 1u"Jy/mas^2",
                    @optic(Tb_peak(_, 5u"GHz")) => 1e11u"K",
                    @optic(visibility(abs, _, SVector(1e8, 0))) => 0.5u"Jy",
                ]
            c = set(c0, o_other, oneunit(o_other(c0)))
            co = modifying(o)(func)
            c_upd = set(c, co, tgt)
            @test coords(c_upd) == coords(c)
            @test o_other(c_upd) == o_other(c)
            @test o(c_upd) != o(c)
            @test func(c_upd) ≈ tgt
        end
        @testset "$o $func" for
                (o, o_other) in [
                    (@optic(_.σ), @optic(_.flux)),
                ],
                (func, tgt) in [
                    fwhm_average => 10u"mas",
                    effective_area => 10u"mas^2",
                ]
            c = set(c0, o_other, oneunit(o_other(c0)))
            co = modifying(o)(func)
            c_upd = set(c, co, tgt)
            @test coords(c_upd) == coords(c)
            @test o_other(c_upd) == o_other(c)
            @test o(c_upd) != o(c)
            @test func(c_upd) ≈ tgt
        end
    end
end

@testitem "uncertainties" begin
    import MonteCarloMeasurements as MCM
    import Uncertain as U
    using StaticArrays
    using Unitful

    uv = SVector(1.0, 0.1)
    uvu = SVector(U.Value(1.0, 0.1), U.Value(0.1, 0.01))
    uvmcm = SVector(MCM.:(±)(1.0, 0.1), MCM.:(±)(0.1, 0.01))

    c = CircularGaussian(flux=U.Value(1.0, 0.1), σ=0.1, coords=SVector(0, 0))
    e = EllipticGaussian(c)
    ec = EllipticGaussianCovmat(c)
    @test fwhm_max(c) == fwhm_min(c) == fwhm_average(c) == 0.23548200450309495
    @test visibility(c, uv) == U.Value(0.8192499856641148 + 0.0im, 0.08192499856641149)
    @test visibility(e, uv) == visibility(c, uv)
    @test visibility(ec, uv) == visibility(c, uv)
    visibility(c, uvu)
    @test_broken (visibility(c, uvmcm); true)

    c = CircularGaussian(flux=U.Value(1.0, 0.1), σ=U.Value(0.1, 0.01), coords=SVector(0, 0))
    e = EllipticGaussian(c)
    ec = EllipticGaussianCovmat(c)
    @test fwhm_max(c) == fwhm_min(c) == fwhm_average(c) == U.Value(0.23548200450309495, 0.023548200450309493)
    @test visibility(c, uv) == U.Value(0.8192499856641148 + 0.0im, 0.08819739670256731)
    @test visibility(e, uv) == U.Value(0.8192499856641148 + 0.0im, 0.08505752517953284)
    @test visibility(ec, uv) == U.Value(0.8192499856641148 + 0.0im, 0.0880787135988441)
    visibility(c, uvu)
    @test_broken (visibility(c, uvmcm); true)

    c = CircularGaussian(flux=MCM.:(±)(1.0, 0.1), σ=0.1, coords=SVector(0, 0))
    e = EllipticGaussian(c)
    ec = EllipticGaussianCovmat(c)
    @test fwhm_max(c) == fwhm_min(c) == fwhm_average(c) == 0.23548200450309495
    # @test visibility(c, uv)
    @test visibility(e, uv) == visibility(c, uv)
    @test visibility(ec, uv) == visibility(c, uv)
    @test_broken (visibility(c, uvu); true)
    visibility(c, uvmcm)
    # @test visibility(c, uvmcm)

    c = CircularGaussian(flux=MCM.:(±)(1.0, 0.1), σ=MCM.:(±)(0.1, 0.01), coords=SVector(0, 0))
    e = EllipticGaussian(c)
    ec = EllipticGaussianCovmat(c)
    @test fwhm_max(c) == fwhm_min(c) == fwhm_average(c) ≈ MCM.:(±)(0.23548200450309495, 0.023548200450309493)
    # @test visibility(c, uv) == MCM.:(±)(0.8192499856641148 + 0.0im, 0.08819739670256731)
    @test visibility(e, uv) ≈ visibility(c, uv)
    @test visibility(ec, uv) ≈ visibility(c, uv)
    @test_broken (visibility(c, uvu); true)
    visibility(c, uvmcm)

    e = EllipticGaussian(flux=U.Value(1.0, 0.1), σ_major=U.Value(0.1, 0.01), ratio_minor_major=U.Value(0.5, 0.05), pa_major=U.Value(0.1, 0.01), coords=SVector(0, 0))
    ec = EllipticGaussianCovmat(e)
    @test visibility(ec, uv) ≈ visibility(e, uv) rtol=1e-3

    e = EllipticGaussian(flux=1.0, σ_major=MCM.:(±)(0.1, 0.01)u"°", ratio_minor_major=MCM.:(±)(0.5, 0.05), pa_major=MCM.:(±)(0.1, 0.01), coords=SVector(0, 0)u"°")
    ec = EllipticGaussianCovmat(e)
    @test visibility(ec, uv) ≈ visibility(e, uv)

    e = EllipticGaussian(flux=MCM.:(±)(1.0, 0.1), σ_major=MCM.:(±)(0.1, 0.01), ratio_minor_major=MCM.:(±)(0.5, 0.05), pa_major=MCM.:(±)(0.1, 0.01), coords=SVector(0, 0))
    ec = EllipticGaussianCovmat(e)
    @test visibility(ec, uv) ≈ visibility(e, uv)

    e = EllipticGaussian(flux=MCM.:(±)(1.0, 0.1), σ_major=MCM.:(±)(0.1, 0.01), ratio_minor_major=MCM.:(±)(0.5, 0.05), pa_major=MCM.:(±)(0.1, 0.01), coords=SVector(MCM.:(±)(0, 0.1), MCM.:(±)(0, 0.1)))
    ec = EllipticGaussianCovmat(e)
    @test visibility(ec, uv) ≈ visibility(e, uv)
end

@testitem "VlbiSkyModels.jl" begin
    import VLBISkyModels as VSM
    using StaticArrays
    using Unitful, UnitfulAngles

    model = VSM.modify(VSM.MRing(0.1, -0.2), VSM.Stretch(VSM.μas2rad(30)))
    @test intensity(model)(SVector(0., 0.)) ≈ 0
    @test intensity(model, SVector(0., 30)u"μas") ≈ 4.514182771546066e20

    @test visibility(model, SVector(0., 0.)) ≈ 1
    @test visibility(model, UV(1, 0.)) ≈ 1
    @test visibility(model, SVector(5e9, 0.)) ≈ -0.304 + 0.0996im  rtol=1e-3
    @test visibility(model)(SVector(5e9, 0.)) ≈ -0.304 + 0.0996im  rtol=1e-3
    @test visibility(abs, model, SVector(5e9, 0.)) ≈ abs(-0.304 + 0.0996im)  rtol=1e-3
end


@testitem "_" begin
    import CompatHelperLocal as CHL
    import Aqua

    CHL.@check()
    Aqua.test_all(InterferometricModels, ambiguities=(;broken=true))
end
