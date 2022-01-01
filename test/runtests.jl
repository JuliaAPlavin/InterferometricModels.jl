using Test
using InterferometricModels
using Unitful, UnitfulAstro
using StaticArrays
import InterferometricModels as IM
import VLBIData as VLBI

# import CompatHelperLocal as CHL
# @testset begin
#     CHL.@check()
# end

@testset begin
    dmod = VLBI.load(VLBI.DifmapModel, "./data/difmap_model.mod")
    mod = IM.model_from_difmap(dmod)
    @test collect(map(flux, components(mod))) ≈ [0.649486, -0.00937331, 1.32054e9]
    @test collect(map(intensity_peak, components(mod))) ≈ [Inf, -0.00022640511888964438, 0.011191592828378834]
    @test collect(map(c -> Tb_peak(c, 15u"GHz"), components(mod))) ≈ [Inf, -1.3934202495658437e6, 6.887914967841405e7]

    @test position_angle(components(mod)[1]) ≈ -0.8165976691060989
    @test_throws MethodError position_angle(components(mod)[2])
    @test position_angle(components(mod)[1] => components(mod)[2]) ≈ 0.510276455129927

    @test_broken visibility(mod, SVector(0, 0))
    @test_broken visibility(mod.components[1], SVector(0, 0))
    @test visibility(mod.components[2], SVector(0, 0)) ≈ -0.00937331 + 0.0im
end
