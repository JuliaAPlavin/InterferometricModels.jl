### A Pluto.jl notebook ###
# v0.17.5

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
	using Revise
	import Pkg
	eval(:(Pkg.develop(path="..")))
	Pkg.resolve()
	using InterferometricModels
end

# ╔═╡ 10d21164-fa86-4af0-a495-8a56932dccba
using VLBIData

# ╔═╡ afc0c2b4-cb50-4eec-8a85-eb8c54ad5d6a
using StaticArrays

# ╔═╡ 50e775c1-625f-4fb0-aab7-31550f7dca1c
using RectiGrids

# ╔═╡ 6980bc53-5338-4e1f-a214-66c49b7602e8
begin
	using PyPlot: plt, matplotlib
	using PyPlotUtils
	pyplot_style!()
end

# ╔═╡ e396f1b3-e22a-4e65-b4ac-c4cdd5886c31
using PlutoUI

# ╔═╡ 1a09ed55-4cf7-43ea-a421-e1d6da45e1c1
using LinearAlgebra

# ╔═╡ 39678544-01fe-43f8-a032-2ba30cc97a43
using PlutoMy

# ╔═╡ 1fb59094-dd9f-479a-a62c-f65a04f40371
using Statistics

# ╔═╡ ee546c1c-39c8-421d-8ff2-cf135eb0cb10
using Difmap

# ╔═╡ 67a0c491-2d7a-4dda-af57-8412fdb6c04a
using BenchmarkTools

# ╔═╡ 2c6f044e-70d7-41e0-a2f4-e4d72fd6c13b
using ProfileSVG

# ╔═╡ 73aa0405-1aab-4321-9c32-9c098c2b8772
using IntervalSets

# ╔═╡ ca9641ba-ca23-4fd2-83a6-b31612a24d09
using Colors

# ╔═╡ 11065b35-d474-4844-af33-4a175881e89e
using Unitful

# ╔═╡ 3c1637e7-44dc-495c-acab-4c0dbb6bd26c
begin
	using Accessors
	ustrip_rad(x::Real) = x  # XXX: may be confusing is passed a unitless model in mas...
	ustrip_rad(x::T) where {T <: Quantity} = dimension(T) == dimension(u"rad") ? ustrip(u"rad", x) : ustrip(x)
	ustrip_rad(x::AbstractArray) = map(ustrip_rad, x)
	ustrip_rad(x::ModelComponent) = @modify(ustrip_rad, x |> Properties())
	ustrip_rad(x::MultiComponentModel) = @modify(ustrip_rad, x.components |> Elements())
end

# ╔═╡ 093ca09d-c09f-48c0-8477-e3cd9c6f9bc9
using DataPipes

# ╔═╡ b61061c6-cbeb-40d4-a09b-5884238a0093


# ╔═╡ fd8d2a74-38a3-405f-81ee-ededc6c71772
mod = VLBI.load("./data/difmap_model.mod")

# ╔═╡ dc475ee0-f5c5-43af-8891-60e4bf834241


# ╔═╡ 821409d1-18da-4116-9c11-3730a1ab5c32
let
	n = 1000
	img = intensity(ustrip(mod)).(grid(SVector, ra=range(10, -10, length=n), dec=range(-10, 10, length=n)))
	plt.figure()
	imshow_ax(img, ColorBar(unit="Jy/mas²"); norm=SymLog(linthresh=1e-3), cmap=:inferno)
	plt.gcf()
end

# ╔═╡ e82945d2-ee22-4e02-98a8-49965f3a73ee
function plot_imageplane(m::MultiComponentModel; kwargs...)
	for c in components(m)
		plot_imageplane(c; kwargs...)
	end
end

# ╔═╡ a0dd1000-ebe1-4fac-83af-475ccc1508a3
function plot_imageplane(c::EllipticGaussian)
	plt.gca().add_patch(
		matplotlib.patches.Ellipse(
			coords(c), width=fwhm_min(c), height=fwhm_max(c), angle=-rad2deg(position_angle(c)),
			fill=false, color=:k, lw=0.5
		)
	)
end

# ╔═╡ 42b78efb-3309-4493-982a-2eaed495e9b2
function plot_imageplane(c::CircularGaussian)
	plt.gca().add_patch(
		matplotlib.patches.Circle(
			coords(c), radius=fwhm_max(c),
			fill=false, color=:k, lw=0.5
		)
	)
end

# ╔═╡ 845dc0fc-d649-4442-aef6-953f0ffb1ecc
let
	plt.figure()
	plot_imageplane(ustrip(mod))
	plt.gca().set_aspect(:equal)
	set_xylims((-5..5)^2; inv=:x)
	xylabels("RA", "Dec"; inline=true)
	plt.gcf()
end

# ╔═╡ 2e00e1fe-1c23-40d1-bf8a-15551d6445db


# ╔═╡ b849ea9e-f703-4d5e-ac43-d0c220869f48


# ╔═╡ 4bb211bf-cf8f-465a-a61c-97194a27454a
mod

# ╔═╡ 8d84c703-bcf7-4d44-897f-bd3834930811
@bind parts MultiCheckBox([:envelopes, :data, :model]; default=[:envelopes, :data, :model])

# ╔═╡ b80db8e3-3753-45d4-a09d-834e1ab1e1f7
mod

# ╔═╡ 94453978-9f72-47d0-9ebd-22d13e8ddbe5
@p VLBI.load("./data/J1033+6051_U_2010_12_24_moj_vis.fits") |> Tables.rows() |> filter(_.stokes ∈ (:RR,)) |> mutate(_.uv.u >= 0 ? (;_.uv, _.visibility) : (;uv=-_.uv, visibility=conj(_.visibility))) |> map(angle(_.visibility) |> rad2deg) |> (length(__), mean(__), std(__), extrema(__))

# ╔═╡ e66c25f5-8f39-474d-86eb-1e518f5760f0
let
	dr = Difmap.execute(
		"""
		observe vis.fits
		select RR
		vis_stats amplitude
		vis_stats phase
		vis_stats real
		vis_stats imaginary
		exit
		""";
		in_files=[
			"./data/J1033+6051_U_2010_12_24_moj_vis.fits" => "vis.fits",
			"./data/difmap_model.mod" => "mod.mod",
		],
	)
	Difmap.inout_pairs(dr)
end

# ╔═╡ cc90ab80-cf2c-4927-9c56-06f683d8f40b
let
	pltnames = ["radplot", "uvplot"]
	tdir = mktempdir()
	dr = Difmap.execute(
		"""
		observe vis.fits
		select RR
		rmodel mod.mod
		rflags="m3"
		device "radplot.ps/VCPS"
		radplot "", 0, 0, 0, 0, -15, +15
		device "uvplot.ps/VCPS"
		uvplot
		exit
		""";
		in_files=[
			"./data/J1033+6051_U_2010_12_24_moj_vis.fits" => "vis.fits",
			"./data/difmap_model.mod" => "mod.mod",
		],
		out_files=["radplot.ps", "uvplot.ps"] .=> tempname.(),
	)
	# Difmap.inout_pairs(dr)
	@assert dr.success
	Difmap.plots(dr)
end

# ╔═╡ d6d4260f-9627-4560-b25a-81b5318d8f57


# ╔═╡ 354ebf8d-d4f4-487d-bb24-d6715e07d81b
data = @p VLBI.load("./data/J1033+6051_U_2010_12_24_moj_vis.fits") |> Tables.rows() |> filter(_.stokes ∈ (:RR,))

# ╔═╡ e50d7eb7-a958-4d90-8a5a-b00bbdc958b7
map(r -> visibility(mod, r.uv), data)

# ╔═╡ 4c25ce58-0138-4ae5-ab5b-2f591ca20a55
begin
	floattype(T::Type, x::Number) = T(x)
	floattype(T::Type, x::AbstractArray) = T.(x)
	floattype(T::Type, x::ModelComponent) = @modify(x->floattype(T,x), x |> Properties())
	floattype(T::Type, x::MultiComponentModel) = @modify(x->floattype(T,x), x.components |> Elements())
end

# ╔═╡ 2345c32d-3db7-42f6-a5ff-a85db528a41c
umod = ustrip_rad(mod)

# ╔═╡ c75ac04d-633a-4313-aa30-ee1c9e68e2ba
map(r -> visibility(umod, r.uv), data)

# ╔═╡ d61f3d23-3c58-4f0b-a46b-96d48b7c168a
map(r -> visibility(umod, r.uv), data)

# ╔═╡ 855e8bc6-c620-4b37-b04f-c23bc0cdfcd5


# ╔═╡ b7cde253-8c16-4a0b-9256-619976f2bf08
@benchmark map(r -> visibility($umod, r.uv), $data)

# ╔═╡ a7cf79e9-e9aa-4daf-8cc0-6e7eb4b64ddc
@benchmark map(r -> $(visibility(mod))(r.uv), $data)

# ╔═╡ 8cdd9b3f-fbb6-4969-91d5-86487b0a5fcb
@benchmark let
	visf = visibility($umod)
	sum(r -> visf(r.uv), $data)
end

# ╔═╡ b41c6a65-7e67-4389-9e2c-20f10504ef91
@profview for _ in 1:1000
	visf = visibility(mod)
	sum(r -> visf(r.uv), data)
end

# ╔═╡ 728e1356-e1f0-408f-9875-e68a52c7f78a


# ╔═╡ 4ef03d45-cb2f-4aa5-9edc-0e9dbeed22dc
begin
	proper_ustrip(::typeof(abs), x) = ustrip(x)
	proper_ustrip(::typeof(angle), x) = ustrip(u"°", x)

	proper_label(::typeof(abs)) = "Amplitude (Jy)"
	proper_label(::typeof(angle)) = "Phase (°)"
end

# ╔═╡ fee5e403-7e44-47e2-bf63-13dea839e029
function radplot(func::typeof(abs), c::ModelComponent, uvrng::AbstractInterval{<:Real}; n=100, alpha=0.2)
	uvdists = range(extrema(uvrng)...; length=n)
	vises = proper_ustrip.(func, visibility_envelope.(func, c, uvdists))
	color = plt.gca()._get_lines.prop_cycler.__next__()["color"] |> mpl_color
	plt.fill_between(uvdists, minimum.(vises), maximum.(vises); ec=color, fc=coloralpha(color, alpha))
	plt.xlabel("UV distance (λ)")
	plt.ylabel(proper_label(func))
	plt.xlim(extrema(uvrng)...)
	func == abs && plt.ylim(0)
end

# ╔═╡ 95883f17-4ccc-4de7-bd29-094e3a8ccb88
function radplot(func::typeof(angle), c::ModelComponent, uvrng::AbstractInterval{<:Real}; n=100, alpha=0.2)
	uvdists = range(extrema(uvrng)...; length=n)
	vises = proper_ustrip.(func, visibility_envelope.(func, c, uvdists))
	color = plt.gca()._get_lines.prop_cycler.__next__()["color"] |> mpl_color
	plt.fill_between(uvdists, minimum.(vises), maximum.(vises); ec=color, fc=coloralpha(color, alpha))
	plt.xlabel("UV distance (λ)")
	plt.ylabel(proper_label(func))

	yl = Interval(plt.ylim()...)
	(0 ± 180) ⊆ yl && plt.ylim(extrema(0 ± 180))	
	plt.xlim(extrema(uvrng)...)
end

# ╔═╡ 541dae20-5a32-452c-b8bf-91854e0aaee1
function radplot(funcs::Vector{<:Function}, args...; ax, kwargs...)
	@assert length(ax) == length(funcs)
	for (i, (a, f)) in zip(ax, funcs) |> enumerate
		plt.sca(a)
		radplot(f, args...; kwargs...)
	end
	for a in ax[begin+1:end-1]
		a.set_xlabel("")
	end
	if length(ax) > 1
		ax[1].xaxis.tick_top()
		ax[1].xaxis.set_label_position(:top)
	end
end

# ╔═╡ 637f8e59-cb27-4929-a127-7beece3059be
begin
	function ylim_shrink(int)
		l = Interval(plt.ylim()...)
		int ⊆ l && plt.ylim(extrema(int))
	end
	
	function ylim_expand(int)
		l = Interval(plt.ylim()...)
		l ⊆ int && plt.ylim(extrema(int))
	end
end

# ╔═╡ 3de9053c-fbf9-4462-ad9c-a45ffca79e6b
function radplot(func, mod, c::Vector{<:NamedTuple})
	uvdists = @p c |> map(norm(_.uv))
	vises = @p c |> map(proper_ustrip(func, func(visibility(mod, _.uv))))
	plt.xlim(auto=true)
	plt.ylim(auto=true)
	plt.scatter(uvdists, vises; s=0.1, c=:red)
	plt.xlabel("UV distance (λ)")
	plt.ylabel(proper_label(func))
	plt.xlim(0)
	func == abs && plt.ylim(0)
	func == angle && ylim_shrink(0±180)
end

# ╔═╡ 12c99ec0-3b30-4bd2-bc75-70075c844a7a
function radplot(func, c::Vector{<:NamedTuple})
	uvdists = @p c |> map(norm(_.uv))
	vises = @p c |> map(proper_ustrip(func, func(_.visibility)))
	plt.xlim(auto=true)
	plt.ylim(auto=true)
	plt.plot(uvdists, vises, ","; color=:k, zorder=1)  # zorder=1 is default for scatter, while zorder=2 for plot
	plt.xlabel("UV distance (λ)")
	plt.ylabel(proper_label(func))
	plt.xlim(0)
	func == abs && plt.ylim(0)
	func == angle && ylim_shrink(0±180)
end

# ╔═╡ 63b0dfb7-3800-44c5-962c-3e005420e993
let mod = VLBI.load("./data/difmap_model.mod")
	data = @p VLBI.load("./data/J1033+6051_U_2010_12_24_moj_vis.fits") |> Tables.rows() |> filter(_.stokes ∈ (:RR,)) |> mutate(_.uv.u >= 0 ? (;_.uv, _.visibility) : (;uv=-_.uv, visibility=conj(_.visibility)))
	_, ax = plt.subplots(2, 1, figsize=(6, 8), gridspec_kw=Dict(:height_ratios => [3, 3], :hspace => 0.1))
	:envelopes ∈ parts && radplot.(Ref([abs, angle]), components(mod), Ref(0..maximum(x->norm(x.uv), data)); ax)
	:data ∈ parts && radplot([abs, angle], data; ax)
	:model ∈ parts && radplot([abs, angle], mod, data; ax)
	plt.ylim(-15, 15)
	plt.gcf()
end

# ╔═╡ 145be160-2b29-4bc1-9a27-4493d787512b


# ╔═╡ b1614a63-ab6f-43de-9b85-25144e3dd4a5


# ╔═╡ f0970c02-75a1-4e14-8173-7d78cdff4fa9


# ╔═╡ 690a88de-95b6-43f5-a616-251f733a00ec
box = matplotlib.offsetbox.AuxTransformBox(plt.gca().transData)
beam = VLBI.image_beam(img)
ellipse = matplotlib.patches.Ellipse((0, 0), ustrip(beam.minor_axis), ustrip(beam.major_axis), -rad2deg(beam.pa), facecolor=:wheat, edgecolor=:k, lw=1) #, hatch="///")
box.add_artist(ellipse)
abox = matplotlib.offsetbox.AnchoredOffsetbox(loc="lower left", pad=0.2, child=box, frameon=false)
plt.gca().add_artist(abox)

# ╔═╡ dd9897e1-e560-4ad0-9396-38cc2a31e8cd


# ╔═╡ a1c927bc-fba3-42f2-9a43-42ffb835da90


# ╔═╡ 10c57c38-ce4a-4753-aeec-5cc2ee85163d
pixelarea(A::KeyedArray) = abs(prod(step.(axiskeys(A))))

# ╔═╡ 0c0b5aa4-5db8-4ade-95fb-49175900ee81
begin
	struct WSlider
	    range::AbstractRange
	    default::Number
	    show_value::Bool
		style::String
	end
	
	
	WSlider(range::AbstractRange; default=missing, show_value=false, width=nothing) = WSlider(range, (default === missing) ? first(range) : default, show_value, isnothing(width) ? nothing : "width: $width")
	
	function Base.show(io::IO, ::MIME"text/html", slider::WSlider)
	    print(io, """<input 
			style="$(slider.style)"
	        type="range" 
	        min="$(first(slider.range))" 
	        step="$(step(slider.range))" 
	        max="$(last(slider.range))" 
	        value="$(slider.default)"
	        $(slider.show_value ? "oninput=\"this.nextElementSibling.value=this.value\"" : "")
	        >""")
	    
	    if slider.show_value
	        print(io, """<output>$(slider.default)</output>""")
	    end
	end
	
	PlutoUI.get(slider::WSlider) = slider.default
end

# ╔═╡ fe1897e2-95a0-4383-8e40-4c529ed15d6f
@bind pa_major° WSlider(-180:5:180; default=0, show_value=true, width="700px")

# ╔═╡ 5bd256aa-02b7-44d2-a45a-2e874f50e110
@bind bsize WSlider(0.1:0.1:5.; show_value=true, width="700px")

# ╔═╡ b475d137-58da-49d5-9fc3-3a91becd29dd
b = beam(EllipticGaussian, σ_major=bsize, ratio_minor_major=0.5, pa_major=deg2rad(15))

# ╔═╡ dcee4f99-f016-4d93-8a91-ffcde40a5259
let
	n = 300
	c = EllipticGaussian(flux=1, σ_major=1, ratio_minor_major=0.5, pa_major=deg2rad(pa_major°), coords=SVector(0, 0))
	figs = []

	int = intensity(convolve(c, b))
	img = int.(grid(SVector, ra=range(10, -10, length=n), dec=range(-10, 10, length=n)))
	plt.figure()
	imshow_ax(img, ColorBar(unit="Jy/b"); norm=SymLog(linthresh=1e-3), cmap=:inferno)
	push!(figs, plt.gcf())

	plt.figure()
	
	slice = int.(grid(SVector, ra=range(-10, 10, length=n), dec=[0]))
	l, = plt.plot(axiskeys(slice, :ra), slice)
	plt.axvline.([-1, 1] .* fwhm_max(c)/2, color=adjust_lightness(l.get_color(), 2))
	
	slice = int.(grid(SVector, ra=[0], dec=range(-10, 10, length=n)))
	l, = plt.plot(axiskeys(slice, :dec), slice')
	plt.axvline.([-1, 1] .* fwhm_min(c)/2, color=adjust_lightness(l.get_color(), 1.5))
	
	push!(figs, plt.gcf())

	figs, maximum(img), sum(img) * pixelarea(img)
end

# ╔═╡ c83145f0-7eb0-446f-bc3e-74f9abf5c041
let
	n = 200
	img = intensity(ustrip(convolve(mod, b))).(grid(SVector, ra=range(10, -10, length=n), dec=range(-10, 10, length=n)))
	plt.figure()
	imshow_ax(img, ColorBar(unit="Jy/b"); norm=SymLog(linthresh=1e-4), cmap=:inferno)
	plt.gcf()
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Accessors = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataPipes = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
Difmap = "697f6f98-b792-4bca-85fe-71b3d053de65"
InterferometricModels = "b395d269-c2ec-4df6-b679-36919ad600ca"
IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoMy = "13de7d54-404c-4927-b51c-7903105debb2"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProfileSVG = "132c30aa-f267-4189-9183-c8a63c7e05e6"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
PyPlotUtils = "5384e752-6c47-47b3-86ac-9d091b110b31"
RectiGrids = "8ac6971d-971d-971d-971d-971d5ab1a71a"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
VLBIData = "679fc9cc-3e84-11e9-251b-cbd013bd8115"

[compat]
Accessors = "~0.1.8"
BenchmarkTools = "~1.2.2"
Colors = "~0.12.8"
DataPipes = "~0.2.5"
Difmap = "~0.1.15"
InterferometricModels = "~0.1.0"
IntervalSets = "~0.5.3"
PlutoMy = "~0.1.9"
PlutoUI = "~0.7.30"
ProfileSVG = "~0.2.1"
PyPlot = "~2.10.0"
PyPlotUtils = "~0.1.3"
RectiGrids = "~0.1.6"
Revise = "~3.3.1"
StaticArrays = "~1.3.2"
Unitful = "~1.10.1"
VLBIData = "~0.3.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.ANSIColoredPrinters]]
git-tree-sha1 = "574baf8110975760d391c710b6341da1afa48d8c"
uuid = "a4c015fc-c6ff-483c-b24f-f7ea428134e9"
version = "0.0.1"

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

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

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
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisKeys]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "IntervalSets", "InvertedIndices", "LazyStack", "LinearAlgebra", "NamedDims", "OffsetArrays", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "d6ff375e8229819ef8e443b25009091e0e88eb24"
uuid = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
version = "0.1.25"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "940001114a0147b6e4d10624276d56d531dd9b49"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.2.2"

[[deps.CFITSIO]]
deps = ["CFITSIO_jll"]
git-tree-sha1 = "4379a2dac795014534b9895a45889aa658fca213"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.4.0"

[[deps.CFITSIO_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg"]
git-tree-sha1 = "2fabb5fc48d185d104ca7ed7444b475705993447"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "3.49.1+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "54fc4400de6e5c3e27be6047da2ef6ba355511f8"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.6"

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
git-tree-sha1 = "6da546248579d3084c06d1babd1431c64f981441"
uuid = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
version = "0.2.5"

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

[[deps.Difmap]]
deps = ["DataPipes", "ImageMagick_jll", "difmap_jll"]
git-tree-sha1 = "f99d8afba929d2685e485c64079483f0484c4912"
uuid = "697f6f98-b792-4bca-85fe-71b3d053de65"
version = "0.1.15"

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

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "67551df041955cc6ee2ed098718c8fcd7fc7aebe"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.12.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.FlameGraphs]]
deps = ["AbstractTrees", "Colors", "FileIO", "FixedPointNumbers", "IndirectArrays", "LeftChildRightSiblingTrees", "Profile"]
git-tree-sha1 = "01c053d887dc6a530a3f1e90b22fc552aae178f2"
uuid = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"
version = "0.2.9"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

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

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f025b79883f361fa1bd80ad132773161d231fd9f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.12+2"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InterferometricModels]]
deps = ["Accessors", "IntervalSets", "LinearAlgebra", "StaticArrays", "Unitful", "UnitfulAstro"]
path = "../../home/aplavin/.julia/dev/InterferometricModels"
uuid = "b395d269-c2ec-4df6-b679-36919ad600ca"
version = "0.1.0"

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
git-tree-sha1 = "22df5b96feef82434b07327e2d3c770a9b21e023"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "3cbe45f4871e60fc142154252322bcf9638c2c1d"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyStack]]
deps = ["LinearAlgebra", "NamedDims", "OffsetArrays", "Test", "ZygoteRules"]
git-tree-sha1 = "a8bf67afad3f1ee59d367267adb7c44ccac7fdee"
uuid = "1fad7336-0346-5a1a-a56f-a06ba010965b"
version = "0.0.7"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "b864cb409e8e445688bc478ef87c0afe4f6d1f8d"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.1.3"

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

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

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
git-tree-sha1 = "f46e8f4e38882b32dcc11c8d31c131d556063f39"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.0"

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

[[deps.MyUnitful]]
deps = ["Unitful", "UnitfulAstro"]
git-tree-sha1 = "9f983d23f225276098e1828e8dd75719a2aafd3f"
uuid = "be63a33b-ca4d-43a5-8045-b0b8c6209429"
version = "0.1.0"

[[deps.NamedDims]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "af6febbfede908c04e19bed954350ac687d892b2"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "0.2.45"

[[deps.Ncurses_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3fb71f3c971fdc8523f9ef3906691124f3048a2"
uuid = "68e3532b-a499-55ff-9963-d1c0c0748b3a"
version = "6.2.0+0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

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

[[deps.PGPLOT_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl", "Pkg"]
git-tree-sha1 = "91bf88f8e171151bdc71e9615159f0a14339f54b"
uuid = "b11e30b1-63be-5002-9df0-88ee0fe906ff"
version = "5.2.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoMy]]
deps = ["ANSIColoredPrinters", "Base64", "BenchmarkTools", "InteractiveUtils", "Markdown", "Suppressor"]
git-tree-sha1 = "cdc0c2863f07665b00dee83b57b9b90c9ee4447b"
uuid = "13de7d54-404c-4927-b51c-7903105debb2"
version = "0.1.9"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5c0eb9099596090bb3215260ceca687b888a1575"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.30"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProfileSVG]]
deps = ["Colors", "FlameGraphs", "Profile", "UUIDs"]
git-tree-sha1 = "e4df82a5dadc26736f106f8d7fc97c42cc6c91ae"
uuid = "132c30aa-f267-4189-9183-c8a63c7e05e6"
version = "0.2.1"

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
deps = ["AxisKeys", "Colors", "DomainSets", "IntervalSets", "OffsetArrays", "PyCall", "PyPlot", "StatsBase", "Unitful"]
git-tree-sha1 = "e1bb1b4ef65e5134121170e14b8cf3f9f232db25"
uuid = "5384e752-6c47-47b3-86ac-9d091b110b31"
version = "0.1.3"

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
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2ae4fe21e97cd13efd857462c1869b73c9f61be3"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.2"

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

[[deps.Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

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
git-tree-sha1 = "dd21b5420bf6e9b76a8c6e56fb575319e7b1f895"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.1"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "c4e1c470a94063b911fd1b1a204cd2bb34a8cd15"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.1.1"

[[deps.VLBIData]]
deps = ["AxisKeys", "DataPipes", "DateFormats", "Dates", "DelimitedFiles", "FITSIO", "InterferometricModels", "MyUnitful", "PyCall", "Reexport", "StaticArrays", "Tables", "Unitful", "UnitfulAstro"]
git-tree-sha1 = "7348d20170c2943a9e71facfa8dda57b82368035"
uuid = "679fc9cc-3e84-11e9-251b-cbd013bd8115"
version = "0.3.1"

[[deps.VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.difmap_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Ncurses_jll", "PGPLOT_jll", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "2c9d61ab9189941639096a9f3dd3396da0d32652"
uuid = "43cd2700-e54a-5645-adf8-b9b288a87fe6"
version = "2.5.11+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═44b64c7e-7621-11ec-33a4-d5b294cc34f6
# ╠═10d21164-fa86-4af0-a495-8a56932dccba
# ╠═afc0c2b4-cb50-4eec-8a85-eb8c54ad5d6a
# ╠═50e775c1-625f-4fb0-aab7-31550f7dca1c
# ╠═6980bc53-5338-4e1f-a214-66c49b7602e8
# ╠═e396f1b3-e22a-4e65-b4ac-c4cdd5886c31
# ╠═b61061c6-cbeb-40d4-a09b-5884238a0093
# ╠═fd8d2a74-38a3-405f-81ee-ededc6c71772
# ╠═dc475ee0-f5c5-43af-8891-60e4bf834241
# ╠═1a09ed55-4cf7-43ea-a421-e1d6da45e1c1
# ╠═39678544-01fe-43f8-a032-2ba30cc97a43
# ╠═fe1897e2-95a0-4383-8e40-4c529ed15d6f
# ╠═5bd256aa-02b7-44d2-a45a-2e874f50e110
# ╠═b475d137-58da-49d5-9fc3-3a91becd29dd
# ╟─dcee4f99-f016-4d93-8a91-ffcde40a5259
# ╟─821409d1-18da-4116-9c11-3730a1ab5c32
# ╟─c83145f0-7eb0-446f-bc3e-74f9abf5c041
# ╠═845dc0fc-d649-4442-aef6-953f0ffb1ecc
# ╠═e82945d2-ee22-4e02-98a8-49965f3a73ee
# ╠═a0dd1000-ebe1-4fac-83af-475ccc1508a3
# ╠═42b78efb-3309-4493-982a-2eaed495e9b2
# ╠═2e00e1fe-1c23-40d1-bf8a-15551d6445db
# ╠═b849ea9e-f703-4d5e-ac43-d0c220869f48
# ╠═4bb211bf-cf8f-465a-a61c-97194a27454a
# ╟─8d84c703-bcf7-4d44-897f-bd3834930811
# ╠═63b0dfb7-3800-44c5-962c-3e005420e993
# ╠═b80db8e3-3753-45d4-a09d-834e1ab1e1f7
# ╠═1fb59094-dd9f-479a-a62c-f65a04f40371
# ╠═94453978-9f72-47d0-9ebd-22d13e8ddbe5
# ╟─e66c25f5-8f39-474d-86eb-1e518f5760f0
# ╟─cc90ab80-cf2c-4927-9c56-06f683d8f40b
# ╠═ee546c1c-39c8-421d-8ff2-cf135eb0cb10
# ╠═d6d4260f-9627-4560-b25a-81b5318d8f57
# ╠═67a0c491-2d7a-4dda-af57-8412fdb6c04a
# ╠═354ebf8d-d4f4-487d-bb24-d6715e07d81b
# ╠═e50d7eb7-a958-4d90-8a5a-b00bbdc958b7
# ╠═3c1637e7-44dc-495c-acab-4c0dbb6bd26c
# ╠═4c25ce58-0138-4ae5-ab5b-2f591ca20a55
# ╠═2345c32d-3db7-42f6-a5ff-a85db528a41c
# ╠═c75ac04d-633a-4313-aa30-ee1c9e68e2ba
# ╠═d61f3d23-3c58-4f0b-a46b-96d48b7c168a
# ╠═855e8bc6-c620-4b37-b04f-c23bc0cdfcd5
# ╠═b7cde253-8c16-4a0b-9256-619976f2bf08
# ╠═a7cf79e9-e9aa-4daf-8cc0-6e7eb4b64ddc
# ╠═8cdd9b3f-fbb6-4969-91d5-86487b0a5fcb
# ╠═2c6f044e-70d7-41e0-a2f4-e4d72fd6c13b
# ╠═b41c6a65-7e67-4389-9e2c-20f10504ef91
# ╠═728e1356-e1f0-408f-9875-e68a52c7f78a
# ╠═73aa0405-1aab-4321-9c32-9c098c2b8772
# ╠═ca9641ba-ca23-4fd2-83a6-b31612a24d09
# ╠═11065b35-d474-4844-af33-4a175881e89e
# ╠═093ca09d-c09f-48c0-8477-e3cd9c6f9bc9
# ╠═4ef03d45-cb2f-4aa5-9edc-0e9dbeed22dc
# ╠═3de9053c-fbf9-4462-ad9c-a45ffca79e6b
# ╠═12c99ec0-3b30-4bd2-bc75-70075c844a7a
# ╠═fee5e403-7e44-47e2-bf63-13dea839e029
# ╠═95883f17-4ccc-4de7-bd29-094e3a8ccb88
# ╠═541dae20-5a32-452c-b8bf-91854e0aaee1
# ╠═637f8e59-cb27-4929-a127-7beece3059be
# ╠═145be160-2b29-4bc1-9a27-4493d787512b
# ╠═b1614a63-ab6f-43de-9b85-25144e3dd4a5
# ╠═f0970c02-75a1-4e14-8173-7d78cdff4fa9
# ╠═690a88de-95b6-43f5-a616-251f733a00ec
# ╠═dd9897e1-e560-4ad0-9396-38cc2a31e8cd
# ╠═a1c927bc-fba3-42f2-9a43-42ffb835da90
# ╠═10c57c38-ce4a-4753-aeec-5cc2ee85163d
# ╠═0c0b5aa4-5db8-4ade-95fb-49175900ee81
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
