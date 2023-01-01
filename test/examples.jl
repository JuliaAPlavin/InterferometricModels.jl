### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ 3c1637e7-44dc-495c-acab-4c0dbb6bd26c
begin
	# used anywhere?..
	using Accessors
	ustrip_rad(x::Real) = x  # XXX: may be confusing is passed a unitless model in mas...
	ustrip_rad(x::T) where {T <: Quantity} = dimension(T) == dimension(u"rad") ? ustrip(u"rad", x) : ustrip(x)
	ustrip_rad(x::AbstractArray) = map(ustrip_rad, x)
	ustrip_rad(x::ModelComponent) = @modify(ustrip_rad, x |> Properties())
	ustrip_rad(x::MultiComponentModel) = @modify(ustrip_rad, x.components |> Elements())
end

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
- Construct a model from a number vector, useful for fitting models
"""

# ╔═╡ 1dc2216d-0103-443a-908a-ee5da9d6cc58


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Accessors = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
DataPipes = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
DisplayAs = "0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"
InterferometricModels = "b395d269-c2ec-4df6-b679-36919ad600ca"
IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PyPlotUtils = "5384e752-6c47-47b3-86ac-9d091b110b31"
RectiGrids = "8ac6971d-971d-971d-971d-971d5ab1a71a"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAngles = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
VLBIData = "679fc9cc-3e84-11e9-251b-cbd013bd8115"
VLBIPlots = "0260e397-8112-41bf-b55a-6b4577718f00"

[compat]
Accessors = "~0.1.8"
DataPipes = "~0.2.7"
DisplayAs = "~0.1.4"
InterferometricModels = "~0.1.1"
IntervalSets = "~0.5.3"
PlutoUI = "~0.7.34"
PyPlotUtils = "~0.1.6"
RectiGrids = "~0.1.6"
Revise = "~3.3.1"
StaticArrays = "~1.3.4"
Unitful = "~1.10.1"
UnitfulAngles = "~0.6.2"
UnitfulAstro = "~1.1.1"
VLBIData = "~0.3.5"
VLBIPlots = "~0.1.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Accessors]]
deps = ["Compat", "CompositionsBase", "ConstructionBase", "Future", "LinearAlgebra", "MacroTools", "Requires", "Test"]
git-tree-sha1 = "0fa53d25794bf1c754f909ee4b0ac31eabff952f"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.8"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1bdcc02836402d104a46f7843b6e6730b1948264"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisKeys]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "IntervalSets", "InvertedIndices", "LazyStack", "LinearAlgebra", "NamedDims", "OffsetArrays", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "d6ff375e8229819ef8e443b25009091e0e88eb24"
uuid = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
version = "0.1.25"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CFITSIO]]
deps = ["CFITSIO_jll"]
git-tree-sha1 = "8425c47db102577eefb93cb37b4480e750116b0d"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.4.1"

[[deps.CFITSIO_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "9c91a9358de42043c3101e3a29e60883345b0b39"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "4.0.0+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f9982ef575e19b0e5c7a98c6e75ee496c0f73a93"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "9aa8a5ebb6b5bf469a7e0e2b5202cf6f8c291104"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.6"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6cdc8832ba11c7695f494c9d9a1c31e90959ce0f"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.6.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.CovarianceEstimation]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "a3e070133acab996660d31dcf479ea42849e368f"
uuid = "587fd27a-f159-11e8-2dae-1979310e6154"
version = "0.2.7"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataPipes]]
deps = ["Accessors", "SplitApplyCombine"]
git-tree-sha1 = "e4b24a30eae6ad440ef909ef48d183a7950d6017"
uuid = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
version = "0.2.7"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DateFormats]]
deps = ["Dates", "DocStringExtensions"]
git-tree-sha1 = "c3f0b713b48dab6e052384e6f4d0e077469bf296"
uuid = "44557152-fe0a-4de1-8405-416d90313ce6"
version = "0.1.12"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dictionaries]]
deps = ["Indexing", "Random"]
git-tree-sha1 = "66bde31636301f4d217a161cabe42536fa754ec8"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.17"

[[deps.DisplayAs]]
git-tree-sha1 = "0cb6c7a4c30a8185cd2a67fdb0d21301bbebbaec"
uuid = "0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"
version = "0.1.4"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.FITSIO]]
deps = ["CFITSIO", "Printf", "Reexport", "Tables"]
git-tree-sha1 = "e6033823834ec0070125120d4d4a1234f1826a47"
uuid = "525bcba6-941b-5504-bd06-fd0dc1a4d2eb"
version = "0.16.12"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InterferometricModels]]
deps = ["Accessors", "IntervalSets", "LinearAlgebra", "StaticArrays", "Unitful", "UnitfulAstro"]
path = "../../home/aplavin/.julia/dev/InterferometricModels"
uuid = "b395d269-c2ec-4df6-b679-36919ad600ca"
version = "0.1.2"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "b55aae9a2bf436fc797d9c253a900913e0e90178"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyStack]]
deps = ["LinearAlgebra", "NamedDims", "OffsetArrays", "Test", "ZygoteRules"]
git-tree-sha1 = "a8bf67afad3f1ee59d367267adb7c44ccac7fdee"
uuid = "1fad7336-0346-5a1a-a56f-a06ba010965b"
version = "0.0.7"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "6b0440822974cab904c8b14d79743565140567f6"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.1"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NamedDims]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "af6febbfede908c04e19bed954350ac687d892b2"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "0.2.45"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.NonNegLeastSquares]]
deps = ["Distributed", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "1271344271ffae97e2855b0287356e6ea5c221cc"
uuid = "b7351bd1-99d9-5c5d-8786-f205a815c4d7"
version = "0.4.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0b5cfbb704034b5b4c1869e36634438a047df065"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8979e9802b4ac3d58c503a20f2824ad67f9074dd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.34"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "71fd4022ecd0c6d20180e23ff1b3e05a143959c2"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.93.0"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "14c1b795b9d764e1784713941e787e1384268103"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.10.0"

[[deps.PyPlotUtils]]
deps = ["Accessors", "AxisKeys", "Colors", "DataPipes", "DomainSets", "IntervalSets", "LinearAlgebra", "NonNegLeastSquares", "OffsetArrays", "PyCall", "PyPlot", "StatsBase", "Unitful"]
git-tree-sha1 = "803be2153093049a12907fae6a44616b101fcbf6"
uuid = "5384e752-6c47-47b3-86ac-9d091b110b31"
version = "0.1.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RectiGrids]]
deps = ["AxisKeys", "Random"]
git-tree-sha1 = "f7e83e8dcc2b6e78d09331720d75cc37355d122f"
uuid = "8ac6971d-971d-971d-971d-971d5ab1a71a"
version = "0.1.6"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "2f9d4d6679b5f0394c52731db3794166f49d5131"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SplitApplyCombine]]
deps = ["Dictionaries", "Indexing"]
git-tree-sha1 = "35efd62f6f8d9142052d9c7a84e35cd1f9d2db29"
uuid = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
version = "1.2.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "d4da8b728580709d736704764e55d6ef38cb7c87"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "a635a9333989a094bddc9f940c04c549cd66afcf"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b95e0b8a8d1b6a6c3e0b3ca393a7a285af47c264"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.10.1"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "d6cfdb6ddeb388af1aea38d2b9905fa014d92d98"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "c4e1c470a94063b911fd1b1a204cd2bb34a8cd15"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.1.1"

[[deps.VLBIData]]
deps = ["AxisKeys", "DataPipes", "DateFormats", "Dates", "DelimitedFiles", "FITSIO", "InterferometricModels", "PyCall", "Reexport", "StaticArrays", "StructArrays", "Tables", "Unitful", "UnitfulAngles", "UnitfulAstro"]
git-tree-sha1 = "e4c331849411fbe1491e86738b31b9d0f06d3e0e"
uuid = "679fc9cc-3e84-11e9-251b-cbd013bd8115"
version = "0.3.5"

[[deps.VLBIPlots]]
deps = ["Colors", "DataPipes", "InterferometricModels", "IntervalSets", "LinearAlgebra", "PyPlotUtils", "Unitful"]
git-tree-sha1 = "f2a51d9775b53ff2378848d7c02ccbbb303f11dd"
uuid = "0260e397-8112-41bf-b55a-6b4577718f00"
version = "0.1.0"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

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
# ╠═3c1637e7-44dc-495c-acab-4c0dbb6bd26c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
