module IntervalArithmeticExt

using InterferometricModels
using InterferometricModels: ustrip, unit
using InterferometricModels: SVector, Interval, endpoints
import IntervalArithmetic as IA


# adapted from RangeEnclosures
function enclose(f, X; maxdepth, atol, rtol)
    fX = f(X)
    if maxdepth == 0 || IA.diam(fX) <= max(atol, rtol * IA.mig(fX))
        return fX
	else
	    Xs = IA.bisect(X)
		Y = mapreduce(IA.hull, Xs) do Xi
			enclose(f, Xi; maxdepth=maxdepth-1, atol, rtol)
		end
	end
end

InterferometricModels.visibility_envelope(::typeof(abs), model, uvdist::Real; maxdepth=10, atol=1e-6, rtol=1e-3) =
	enclose(
		θ -> visibility(abs, model, uvdist * SVector(sincos(θ))) |> ustrip,
		IA.interval(0, π);
		maxdepth, atol, rtol
	) |> _Interval |> i->_mul(i, unit(flux(model)))

InterferometricModels.visibility_envelope(::typeof(angle), model, uvdist::Real; maxdepth=10, atol=1e-6, rtol=1e-3) =
	enclose(
		θ -> visibility(angle, model, uvdist * SVector(sincos(θ))),
		IA.interval(0, 2π);
		maxdepth, atol, rtol
	) |> _Interval


# XXX: should be an extension in IntervalArithmetic
_Interval(i::IA.Interval) = Interval(IA.inf(i), IA.sup(i))

# XXX
_mul(x::Interval, u) = Interval((endpoints(ustrip(x)) .* u)...)

end
