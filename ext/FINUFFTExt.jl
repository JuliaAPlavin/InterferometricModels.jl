module FINUFFTExt

import FINUFFT
import InterferometricModels
using InterferometricModels: MultiComponentModel, Point, UVType, components, coords, flux, visibility
using Unitful: NoUnits, unit, ustrip

# Σⱼ Fⱼ·cis(2π uv·xⱼ) at every uv point as one type-3 NUFFT, to relative tolerance 1e-9
function InterferometricModels.visibilities(m::MultiComponentModel{<:AbstractVector{<:Point}}, uvs::AbstractVector{<:UVType})
    # FINUFFT corrupts memory when given no target points
    isempty(uvs) && return visibility.(Ref(m), uvs)
    cs = components(m)
    # units go into the coordinates, so that uv·x is plain; FINUFFT takes Float64 arrays only
    uvunit = oneunit(eltype(eltype(uvs)))
    x, y = (map(c -> Float64(2π * NoUnits(coords(c)[k] * uvunit)), cs) for k in 1:2)
    F = flux.(cs)
    V = FINUFFT.nufft2d3(x, y, complex.(Float64.(ustrip.(F))), +1, 1e-9,
                         map(uv -> Float64(ustrip(uv[1])), uvs), map(uv -> Float64(ustrip(uv[2])), uvs);
                         nthreads = Threads.nthreads())
    V .* unit(eltype(F))
end

end
