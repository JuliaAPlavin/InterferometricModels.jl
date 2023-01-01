module Transformations

abstract type Transformation end
function transform end

struct Identity <: Transformation end
transform(tr::Identity, x, ::Val) = x

struct FuncTransformation{T1, T2, T3} <: Transformation
	flux::T1
	σ::T2
	coords::T3
end

function FuncTransformation(; flux, σ, x=nothing, coords=nothing)
    @assert isnothing(x) != isnothing(coords)
    coords = @something(coords, c -> (x(c[1]), x(c[2])))
    return FuncTransformation(flux, σ, coords)
end

transform(tr::FuncTransformation, x, ::Val{:flux}) = tr.flux(x)
transform(tr::FuncTransformation, x, ::Val{:σ}) = tr.σ(x)
transform(tr::FuncTransformation, x, ::Val{:coords}) = tr.coords(x)
end


n_params(::Type{<:CircularGaussian}) = 4
function from_vec(::Type{<:CircularGaussian}, x::AbstractVector, start=1)
    CircularGaussian(x[start], x[start+1], SVector(x[start+2], x[start+3]))
end
function transform_vec!(::Type{<:CircularGaussian}, tr::Transformations.Transformation, x::AbstractVector, start=1)
    x[start] = Transformations.transform(tr, x[start], Val(:flux))
    x[start+1] = Transformations.transform(tr, x[start+1], Val(:σ))
    x[start+2], x[start+3] = Transformations.transform(tr, SVector(x[start+2], x[start+3]), Val(:coords))
end
param_names(::Type{<:CircularGaussian}) = ["flux", "σ", "x", "y"]


@generated n_params(::Type{MultiComponentModel{TUP}}) where {TUP <: Tuple} = sum(n_params(T) for T in fieldtypes(TUP))
@generated function from_vec(::Type{MultiComponentModel{TUP}}, x::AbstractVector, start=1) where {TUP <: Tuple}
    types = fieldtypes(TUP)
    quote
        @assert length(x) == $(n_params(MultiComponentModel{TUP}))
        MultiComponentModel(tuple(
            $([:(from_vec($T, x, start + $(sum(n_params(TT) for TT in types[1:i-1]; init=0)))) for (i, T) in enumerate(types)]...)
        ))
    end
end

@generated function transform_vec!(::Type{MultiComponentModel{TUP}}, tr::Transformations.Transformation, x::AbstractVector, start=1) where {TUP <: Tuple}
    types = fieldtypes(TUP)
    quote
        @assert length(x) == $(n_params(MultiComponentModel{TUP}))
        $([:(transform_vec!($T, tr, x, start + $(sum(n_params(TT) for TT in types[1:i-1]; init=0)))) for (i, T) in enumerate(types)]...)
        return x
    end
end
transform_vec(TMOD::Type{<:MultiComponentModel}, tr::Transformations.Transformation, x::AbstractVector) = transform_vec!(TMOD, tr, copy(x))

param_names(::Type{MultiComponentModel{TUP}}) where {TUP <: Tuple} = mapreduce(vcat, fieldtypes(TUP) |> enumerate) do (i, T)
    ["$(i)_$(n)" for n in param_names(T)]
end
