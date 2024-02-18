using AxisKeys
using NFFT, FINUFFT
using DataPipes

include(joinpath(dirname(pathof(NFFT)), "..", "Wrappers", "FINUFFT.jl"))

export ImageVisibilityOperator

struct ImageVisibilityOperator{IM<:KeyedArray{<:Complex},VI<:Vector{<:Complex},UV,P,AP}
    img₀::IM
    vis₀::VI
    vis_muls::VI
    uv::UV
    plan::P
    adj_plan::AP
end

Base.adjoint(vo::ImageVisibilityOperator) = Adjoint(vo)

function ImageVisibilityOperator(img::KeyedArray, uv::Vector{<:UVType}; reltol=1e-10)
    steps = @p axiskeys(img) .|> step
    @assert allequal(map(abs, steps))
    uvmat = @p let
        uv
        map(_ .* ustrip.(u"rad", steps))
        stack
        -  # XXX: why negate? equivalent to complex conjugate of visibilities
    end

    offset = map(axiskeys(img)) do aks
        @assert iseven(length(aks))
        ustrip(u"rad", aks[length(aks) ÷ 2 + 1])
    end
    vis_muls = cispi.(2 .* dot.(Ref(offset), uv))

    p = FINUFFTPlan(uvmat, size(img); reltol)
    ImageVisibilityOperator(copy(img), zeros(eltype(img), length(uv)), vis_muls, uv, p, adjoint(p))
end

function Base.:*(p::ImageVisibilityOperator, img::KeyedArray{<:Complex}; kargs...)
    vis = similar(img, length(p.uv))
    return mul!(vis, p, img; kargs...)
end

function Base.:*(p::Adjoint{<:Any,<:ImageVisibilityOperator}, vis::AbstractVector; kargs...)
    return mul!(copy(p.parent.img₀), p, vis; kargs...)
end

function LinearAlgebra.mul!(vis, p::ImageVisibilityOperator, img::KeyedArray; kargs...)
    mul!(vis, p.plan, AxisKeys.keyless_unname(img); kargs...)
    vis .*= p.vis_muls
    return vis
end

function LinearAlgebra.mul!(img, p::Adjoint{<:Any,<:ImageVisibilityOperator}, vis; kargs...)
    p.parent.vis₀ .= vis .* conj.(p.parent.vis_muls)
    mul!(AxisKeys.keyless_unname(img), p.parent.adj_plan, p.parent.vis₀; kargs...)
    return img
end
