using Test
using InterferometricModels
using Unitful, UnitfulAstro
using Utils
import InterferometricModels as IM
import VLBIData as VLBI

# import CompatHelperLocal as CHL
# @testset begin
#     CHL.@check()
# end

@testset begin
    dmod = VLBI.load(VLBI.DifmapModel, "./data/difmap_model.mod")
    mod = IM.model_from_difmap(dmod)
    @test collect(map(flux, components(mod))) ≈ [0.649486u"Jy", -0.00937331u"Jy", 1.32054e9u"Jy"]
    @test collect(map(intensity_peak, components(mod))) ≈ [Inf*u"Jy*mas^-2", -0.00022640511888964438u"Jy*mas^-2", 0.011191592828378834u"Jy*mas^-2"]
    @test collect(map(c -> Tb_peak(c, 15u"GHz"), components(mod))) ≈ [Inf*u"K", -1.3934202495658437e6u"K", 6.887914967841405e7u"K"]
end
