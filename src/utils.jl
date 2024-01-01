σ_to_fwhm(σ) = √(8 * log(2)) * σ
fwhm_to_σ(fwhm) = fwhm / √(8 * log(2))
Accessors.inverse(::typeof(σ_to_fwhm)) = fwhm_to_σ
Accessors.inverse(::typeof(fwhm_to_σ)) = σ_to_fwhm


# unitful: use units
intensity_to_Tb(intensity, ν) = intensity * u"c"^2 / (2 * u"k" * ν^2) |> u"K"
# plain number: assume janskys/mas^2
intensity_to_Tb(intensity::Real, ν) = ustrip(u"K", intensity_to_Tb(intensity*u"Jy/(1e-3*arcsecond)^2", ν))
