### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 44b64c7e-7621-11ec-33a4-d5b294cc34f6
begin
	# using Revise
	# import Pkg
	# eval(:(Pkg.develop(path="..")))
	# Pkg.resolve()
	using InterferometricModels
end

# ╔═╡ 208d6741-dc9b-4ccf-ac7c-58b0fa8a9999
using LinearAlgebra: norm

# ╔═╡ a1fb64b7-cc8a-463c-b457-b5dac2e2aed4
using PyPlotUtils; pyplot_style!()

# ╔═╡ c3debbf2-eff4-4284-87d6-cfc2be70bfbc
using VLBIPlots

# ╔═╡ 10d21164-fa86-4af0-a495-8a56932dccba
using VLBIData

# ╔═╡ afc0c2b4-cb50-4eec-8a85-eb8c54ad5d6a
using StaticArrays

# ╔═╡ 50e775c1-625f-4fb0-aab7-31550f7dca1c
using RectiGrids

# ╔═╡ e396f1b3-e22a-4e65-b4ac-c4cdd5886c31
using PlutoUI

# ╔═╡ 19fa40a7-581d-4b75-bcfe-82a1b8fc6c02
using DisplayAs: Text as AsText

# ╔═╡ 73aa0405-1aab-4321-9c32-9c098c2b8772
using IntervalSets

# ╔═╡ 11065b35-d474-4844-af33-4a175881e89e
using Unitful, UnitfulAstro, UnitfulAngles

# ╔═╡ 093ca09d-c09f-48c0-8477-e3cd9c6f9bc9
using DataPipes

# ╔═╡ 606df552-7a32-41a7-ae87-0c381adc85d5
using Accessors

# ╔═╡ b61061c6-cbeb-40d4-a09b-5884238a0093
md"""
!!! info "InterferometricModels.jl"
	`InterferometricModels` defines model components typically used to describe structure visible in VLBI observations.

These model components currently include point sources (`Point`), circular and elliptic gaussian functions (`CircularGaussian`, `EllipticGaussian`). Components can be combined into a composite `MultiComponentModel`.

Components and models support computing their intensity (image plane) and interferometric visibility (UV plane). Some other, more advanced or convenient functions, are illustrated in examples below.
"""

# ╔═╡ 1353a730-e416-40c8-b86d-72c7503ef375
md"""
# Simple component
"""

# ╔═╡ dc475ee0-f5c5-43af-8891-60e4bf834241
md"""
First, let's create a simple single-component model (dimensionless for now) with its position angle controlled by a slider:
"""

# ╔═╡ fe1897e2-95a0-4383-8e40-4c529ed15d6f
md"""
PA: $(@bind pa_major Slider((-180:5:180)u"°"; default=0, show_value=true))
"""

# ╔═╡ aae6e366-7856-4560-8c90-eb7e9a7d0ec2
comp = EllipticGaussian(flux=1, σ_major=0.5, ratio_minor_major=0.5, pa_major=convert(Float64, pa_major), coords=SVector(0., 0))

# ╔═╡ 83bce173-f6f1-4a56-9fd4-caee52fa65fe
md"""
... and compute its intensity image:
"""

# ╔═╡ 46c40ba3-a0a9-4385-8d8f-1b0630574d87
comp_img = let
	n = 200
	intensity(comp).(grid(SVector, ra=range(10, -10, length=n), dec=range(-10, 10, length=n)))
end;

# ╔═╡ 3bc23b30-ee51-4591-b803-42d9c07cf978
comp_img |> AsText

# ╔═╡ 9c7cb933-bef4-444a-a502-b68969a3fa47
let
	plt.figure()
	imshow_ax(comp_img, ColorBar(); norm=SymLog(linthresh=1e-3), cmap=:inferno)
	xylabels(comp_img; inline=true)
	plt.gcf()
end

# ╔═╡ c1c73904-777d-42d3-a27c-d92c38b18ff2
md"""
The intensity peak value is often useful, and can be computed directly:
"""

# ╔═╡ 767d69fb-b30b-4c14-a07c-819818b447c4
intensity_peak(comp), maximum(comp_img)

# ╔═╡ c21afcce-6631-485d-bc8c-4cf11bdd7489
md"""
One-dimensional intensity slices can directly demonstrate the meaning of `FWHM` (full width at half maximum). In the following plot, max and min FWHM values are marked. They should agree with the computed curve depending on the positional angle selected above:
"""

# ╔═╡ 11b885ac-cb5a-4f0f-9ee0-997a85a6d037
let
	n = 200

	plt.figure()
	
	slice = intensity(comp).(grid(SVector, ra=range(-2, 2, length=n), dec=[0]))
	l, = plt.plot(axiskeys(slice, :ra), slice, color=:C1)
	plt.axvline.([-1, 1] .* fwhm_max(comp)/2, color=adjust_lightness(:C0, 1.5))
	plt.axvline.([-1, 1] .* fwhm_min(comp)/2, color=adjust_lightness(:C2, 1.5))
	
	plt.gcf()
end

# ╔═╡ 59a48752-f17f-4ed0-a061-e1474f3e310c
md"""
We can also convolve the component with another one - typically, a "beam". Here, the beam is also an elliptical gaussian function with a different orientation:
"""

# ╔═╡ 5bd256aa-02b7-44d2-a45a-2e874f50e110
md"""
Beam size: $(@bind bsize Slider(0.1:0.1:5.; show_value=true, default=2))
"""

# ╔═╡ b475d137-58da-49d5-9fc3-3a91becd29dd
b = beam(EllipticGaussian, σ_major=bsize, ratio_minor_major=0.5, pa_major=deg2rad(15))

# ╔═╡ f32b584a-07d9-4e47-b694-b7c5fc6842ee
md"""
`beam()` is a simple convenience constructor: it sets flux so that the component is properly normalized to have the intensity peak of unity, and has zero coordinates.
"""

# ╔═╡ 8bc66f52-5c8d-4c52-a310-9db4e8d5259b
comp_convolved = convolve(comp, b)

# ╔═╡ dcee4f99-f016-4d93-8a91-ffcde40a5259
let
	n = 200
	img = intensity(comp_convolved).(grid(SVector, ra=range(10, -10, length=n), dec=range(-10, 10, length=n)))
	plt.figure()
	imshow_ax(img, ColorBar(unit="/b"); norm=SymLog(linthresh=1e-3), cmap=:inferno)
	xylabels(img; inline=true)
	plt.gcf()
end

# ╔═╡ cb89545d-94f2-4268-8ce1-43a4629a883f
md"""
# Real multicomponent model
"""

# ╔═╡ 3072d632-197d-495b-9fd9-3e6a5e1ffcec
md"""
Now let's read a real source model from a file, using the `VLBIData.jl` package for IO:
"""

# ╔═╡ fd8d2a74-38a3-405f-81ee-ededc6c71772
mod = VLBI.load("./data/difmap_model.mod")

# ╔═╡ 3cf4f288-f78f-4a14-8d8a-0d99ed35e737
md"""
You can inspect each model component and its parameters interactively above.

Note that this model has proper units: flux in `Janskys`, coordinates and sizes in `mas`. Sometimes raw float numbers are needed, and `InterferometricModels` defines corresponding `ustrip` methods:
"""

# ╔═╡ a0d35df3-bd30-468e-962a-ba15537db935
ustrip(mod)

# ╔═╡ d8bb3949-ee87-4c23-9276-400335d1f113
md"""
This is a simple visualization of all three model components, showing their locations and sizes:
"""

# ╔═╡ 845dc0fc-d649-4442-aef6-953f0ffb1ecc
let
	plt.figure()
	plot_imageplane(ustrip(mod))
	plt.gca().set_aspect(:equal)
	set_xylims((-5..5)^2; inv=:x)
	xylabels("RA", "Dec"; inline=true)
	plt.gcf()
end

# ╔═╡ cfd9cd8e-fc72-4f76-824b-8fef606fab5e
md"""
Compute the intensity image of the whole model:
"""

# ╔═╡ 8798c312-5fad-437d-9c1b-0139bff875f2
mod_img = let
	n = 1000
	img = intensity(mod).(grid(SVector, ra=range(10, -10, length=n)u"mas", dec=range(-10, 10, length=n)u"mas"))
end;

# ╔═╡ 5ddb5d60-100a-46d6-a738-474382563d50
mod_img |> AsText

# ╔═╡ 4ad7df4c-7604-499c-9045-01e7e59a4ba0
md"""
It also has proper units, both in values and in axes.
"""

# ╔═╡ c3eafe4b-9183-4bc7-9cba-7325ed530af2
md"""
Here is plot of the intensity image:
"""

# ╔═╡ 821409d1-18da-4116-9c11-3730a1ab5c32
let
	plt.figure()
	imshow_ax(ustrip.(mod_img), ColorBar(unit="Jy/mas²"); norm=SymLog(linthresh=1e-3), cmap=:inferno)
	xylabels(mod_img; inline=true)
	plt.gcf()
end

# ╔═╡ f2c4d92d-6979-484a-b289-83b005152124
md"""
And the intensity image convolved with the beam we defined above:
"""

# ╔═╡ c83145f0-7eb0-446f-bc3e-74f9abf5c041
let
	n = 200
	img = intensity(ustrip(convolve(mod, b))).(grid(SVector, ra=range(10, -10, length=n), dec=range(-10, 10, length=n)))
	plt.figure()
	imshow_ax(img, ColorBar(unit="Jy/b"); norm=SymLog(linthresh=1e-4), cmap=:inferno)
	xylabels("RA", "Dec"; inline=true)
	plt.gcf()
end

# ╔═╡ 3a55425c-f6e0-45a3-97da-6ff9c611682c
md"""
`InterferometricModels` has features specific to interferometry: most notably, computation of the visibility function.

First, let's compute the model visibility at a few points.
"""

# ╔═╡ e9888326-9b3f-4ed6-b714-f64e4349c87a
md"""
At zero, should equal to the total flux in the model:
"""

# ╔═╡ f5cbfafb-d495-4ccb-be80-6774f4571aed
visibility(mod)(SVector(0, 0)), flux(mod)

# ╔═╡ b7373f28-bc61-4f27-87d7-004caa75620f
md"""
At a very large UV distance, should be close to zero:
"""

# ╔═╡ fd0f4511-61e4-44dc-b0c3-923c32898ecc
visibility(mod)(SVector(3e10, 0))

# ╔═╡ 26e114f7-0163-4b89-816f-7dd59e5accd8
md"""
At a reasonable UV distance, corresponding to the model scale:
"""

# ╔═╡ 94cbdc30-9715-4502-b8b6-7538e1894172
visibility(mod)(SVector(2e9, 0))

# ╔═╡ f34a232c-4946-41fd-af30-87a6e0e5cee4
md"""
And at a set of random points on the UV plane:
"""

# ╔═╡ 63487efd-72ca-492b-babb-5c9fef2619c9
uvs = [1e9 * SVector(randn(2)...) for _ in 1:1000]

# ╔═╡ bc5a31db-3619-4c17-80f0-7e561eebb7b4
vises = visibility(mod).(uvs)

# ╔═╡ 580e9844-aa32-48a2-b949-46d695f373a1
let
	plt.figure()
	plt.scatter(norm.(uvs), abs.(ustrip.(vises)); s=1, color=:k)
	xylabels("UV distance (λ)", "Visibility amplitude (Jy)")
	plt.gcf()
end

# ╔═╡ 5dbfa878-ebd2-4112-bc98-0f2921c2c7f5
md"""
It's also convenient to compare contributions of each model component at different UV distances. `InterferometricModels` provides the `visibility_envelope` function: it calculates the range of possible visibility amplitudes and phases of a single component at a fixed UV distance:
"""

# ╔═╡ e6929ddf-3236-4627-be5e-d5db9c2a3bc0
visibility_envelope(abs, components(mod)[1], 1e9)

# ╔═╡ 4d2a6357-889a-4ef9-b1b0-739693f5cf33
md"""
Let's plot envelopes for all three components in our model:
"""

# ╔═╡ 49be2200-42f7-4437-ad8e-b4894c501652
let
	_, ax = plt.subplots(2, 1, gridspec_kw=Dict(:height_ratios => [3, 1]))
	for (i, c) in components(mod) |> enumerate
		radplot([abs, angle], c, 0..1e9; ax, label="Comp #$i")
	end
	plt.sca(ax[1])
	legend_inline_right()
	plt.gcf()
end

# ╔═╡ 34d8af3e-f43f-480d-9d1d-1d22372202ac
md"""
# Neat example
"""

# ╔═╡ 28054c1d-4612-4dfc-bf8b-73b3709caa9b
md"""
Let's create a simple model and visualize how its complex visibility function behaves with varying UV location.
"""

# ╔═╡ eeac5754-92a8-40a3-a483-4ddc9dcbb4ce
emod = MultiComponentModel((
	CircularGaussian(flux=1u"Jy", σ=0.1u"mas", coords=SVector(0, 0.03)u"mas"),
	CircularGaussian(flux=0.7u"Jy", σ=0.3u"mas", coords=SVector(0, 2.)u"mas"),
	CircularGaussian(flux=0.2u"Jy", σ=0.7u"mas", coords=SVector(0, 5.)u"mas"),
))

# ╔═╡ 3ae7bb8d-e8e0-4b3a-a625-04ab72aba2f4
md"""
UV location is controlled with two slides, setting distance and angle:
"""

# ╔═╡ be50bb46-7110-4ae5-932b-bd5fd7601010
md"""
UV dist: $(@bind uvdist Slider(0:2e6:1e8; default=1e7, show_value=true))

UV angle: $(@bind uvang Slider((-90:5:90)u"°"; default=30u"°", show_value=true))
"""

# ╔═╡ 3fec21ab-9497-4603-bd61-0ca398e5476d
uv = uvdist .* SVector(sincos(uvang)) |> reverse

# ╔═╡ 3c2254aa-9fc8-489a-b368-897cd2cb603b
md"""
And plot the complex visibility function at this specific UV location:
"""

# ╔═╡ 3c5d095a-d767-43d9-9fe9-da3183d079eb
let
	plt.figure()
	plt.gca().set_aspect(:equal)

	plt.gca().add_artist(matplotlib.patches.Circle(xy=(0, 0), radius=ustrip(flux(emod)), fc=:none, ec=:grey))
	plt.text(1.7, -0.8, "Total flux", ha=:right, color=:grey)

	map(0:length(components(emod))) do n
		curvis = @p sum(visibility(components(emod)[i], uv) for i in 1:n; init=0u"Jy") |> ustrip(u"Jy", __)
		prevvis = @p sum(visibility(components(emod)[i], uv) for i in 1:(n-1); init=0u"Jy") |> ustrip(u"Jy", __)
		n > 0 && plt.plot(real.([prevvis, curvis]), imag.([prevvis, curvis]), label=n == 0 ? "Total" : "Comp #$n")
		n == length(components(emod)) && plt.plot([0, real(curvis)], [0, imag(curvis)]; color=:k, ls=":")
		plt.scatter(real(curvis), imag(curvis); color=n == 0 ? :k : nothing, s=2, zorder=10, label=n == 0 ? "Total" : nothing)
	end

	plt.legend()
	set_xylims((-0.1..2) × (0±1))
	xylabels("Re (Jy)", "Im (Jy)"; inline=false)
	plt.gcf()
end

# ╔═╡ d6d4260f-9627-4560-b25a-81b5318d8f57
md"""
# Advanced functionality
"""

# ╔═╡ 43d32d88-a626-4ebf-b232-5c4533d63554
md"""
Currently undocumented, see source or tests:
- Compute brightness temperatures
- `Accessors.jl` integration
- and a bit more
"""

# ╔═╡ 1dc2216d-0103-443a-908a-ee5da9d6cc58


# ╔═╡ 3c1637e7-44dc-495c-acab-4c0dbb6bd26c
begin
	# used anywhere?..
	ustrip_rad(x::Real) = x  # XXX: may be confusing is passed a unitless model in mas...
	ustrip_rad(x::T) where {T <: Quantity} = dimension(T) == dimension(u"rad") ? ustrip(u"rad", x) : ustrip(x)
	ustrip_rad(x::AbstractArray) = map(ustrip_rad, x)
	ustrip_rad(x::ModelComponent) = @modify(ustrip_rad, x |> Properties())
	ustrip_rad(x::MultiComponentModel) = @modify(ustrip_rad, x.components |> Elements())
end

# ╔═╡ Cell order:
# ╟─b61061c6-cbeb-40d4-a09b-5884238a0093
# ╟─1353a730-e416-40c8-b86d-72c7503ef375
# ╟─dc475ee0-f5c5-43af-8891-60e4bf834241
# ╟─fe1897e2-95a0-4383-8e40-4c529ed15d6f
# ╠═aae6e366-7856-4560-8c90-eb7e9a7d0ec2
# ╟─83bce173-f6f1-4a56-9fd4-caee52fa65fe
# ╟─3bc23b30-ee51-4591-b803-42d9c07cf978
# ╠═46c40ba3-a0a9-4385-8d8f-1b0630574d87
# ╠═9c7cb933-bef4-444a-a502-b68969a3fa47
# ╟─c1c73904-777d-42d3-a27c-d92c38b18ff2
# ╠═767d69fb-b30b-4c14-a07c-819818b447c4
# ╟─c21afcce-6631-485d-bc8c-4cf11bdd7489
# ╠═11b885ac-cb5a-4f0f-9ee0-997a85a6d037
# ╟─59a48752-f17f-4ed0-a061-e1474f3e310c
# ╟─5bd256aa-02b7-44d2-a45a-2e874f50e110
# ╠═b475d137-58da-49d5-9fc3-3a91becd29dd
# ╟─f32b584a-07d9-4e47-b694-b7c5fc6842ee
# ╠═8bc66f52-5c8d-4c52-a310-9db4e8d5259b
# ╟─dcee4f99-f016-4d93-8a91-ffcde40a5259
# ╟─cb89545d-94f2-4268-8ce1-43a4629a883f
# ╟─3072d632-197d-495b-9fd9-3e6a5e1ffcec
# ╠═fd8d2a74-38a3-405f-81ee-ededc6c71772
# ╟─3cf4f288-f78f-4a14-8d8a-0d99ed35e737
# ╠═a0d35df3-bd30-468e-962a-ba15537db935
# ╟─d8bb3949-ee87-4c23-9276-400335d1f113
# ╟─845dc0fc-d649-4442-aef6-953f0ffb1ecc
# ╟─cfd9cd8e-fc72-4f76-824b-8fef606fab5e
# ╟─5ddb5d60-100a-46d6-a738-474382563d50
# ╠═8798c312-5fad-437d-9c1b-0139bff875f2
# ╟─4ad7df4c-7604-499c-9045-01e7e59a4ba0
# ╟─c3eafe4b-9183-4bc7-9cba-7325ed530af2
# ╠═821409d1-18da-4116-9c11-3730a1ab5c32
# ╟─f2c4d92d-6979-484a-b289-83b005152124
# ╟─c83145f0-7eb0-446f-bc3e-74f9abf5c041
# ╟─3a55425c-f6e0-45a3-97da-6ff9c611682c
# ╟─e9888326-9b3f-4ed6-b714-f64e4349c87a
# ╠═f5cbfafb-d495-4ccb-be80-6774f4571aed
# ╟─b7373f28-bc61-4f27-87d7-004caa75620f
# ╠═fd0f4511-61e4-44dc-b0c3-923c32898ecc
# ╟─26e114f7-0163-4b89-816f-7dd59e5accd8
# ╠═94cbdc30-9715-4502-b8b6-7538e1894172
# ╟─f34a232c-4946-41fd-af30-87a6e0e5cee4
# ╠═63487efd-72ca-492b-babb-5c9fef2619c9
# ╠═bc5a31db-3619-4c17-80f0-7e561eebb7b4
# ╟─580e9844-aa32-48a2-b949-46d695f373a1
# ╟─5dbfa878-ebd2-4112-bc98-0f2921c2c7f5
# ╠═e6929ddf-3236-4627-be5e-d5db9c2a3bc0
# ╟─4d2a6357-889a-4ef9-b1b0-739693f5cf33
# ╠═49be2200-42f7-4437-ad8e-b4894c501652
# ╟─34d8af3e-f43f-480d-9d1d-1d22372202ac
# ╟─28054c1d-4612-4dfc-bf8b-73b3709caa9b
# ╠═eeac5754-92a8-40a3-a483-4ddc9dcbb4ce
# ╟─3ae7bb8d-e8e0-4b3a-a625-04ab72aba2f4
# ╟─be50bb46-7110-4ae5-932b-bd5fd7601010
# ╠═3fec21ab-9497-4603-bd61-0ca398e5476d
# ╟─3c2254aa-9fc8-489a-b368-897cd2cb603b
# ╟─3c5d095a-d767-43d9-9fe9-da3183d079eb
# ╟─d6d4260f-9627-4560-b25a-81b5318d8f57
# ╟─43d32d88-a626-4ebf-b232-5c4533d63554
# ╠═1dc2216d-0103-443a-908a-ee5da9d6cc58
# ╠═44b64c7e-7621-11ec-33a4-d5b294cc34f6
# ╠═208d6741-dc9b-4ccf-ac7c-58b0fa8a9999
# ╠═a1fb64b7-cc8a-463c-b457-b5dac2e2aed4
# ╠═c3debbf2-eff4-4284-87d6-cfc2be70bfbc
# ╠═10d21164-fa86-4af0-a495-8a56932dccba
# ╠═afc0c2b4-cb50-4eec-8a85-eb8c54ad5d6a
# ╠═50e775c1-625f-4fb0-aab7-31550f7dca1c
# ╠═e396f1b3-e22a-4e65-b4ac-c4cdd5886c31
# ╠═19fa40a7-581d-4b75-bcfe-82a1b8fc6c02
# ╠═73aa0405-1aab-4321-9c32-9c098c2b8772
# ╠═11065b35-d474-4844-af33-4a175881e89e
# ╠═093ca09d-c09f-48c0-8477-e3cd9c6f9bc9
# ╠═606df552-7a32-41a7-ae87-0c381adc85d5
# ╠═3c1637e7-44dc-495c-acab-4c0dbb6bd26c
