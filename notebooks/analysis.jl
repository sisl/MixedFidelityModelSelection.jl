### A Pluto.jl notebook ###
# v0.17.7

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

# ╔═╡ 9bb8e120-8469-11ec-11ba-6d196236bb2f
using BSON

# ╔═╡ fb741838-9f6d-4920-bc59-c89e01d20d99
using Plots; default(fontfamily="Computer Modern", framestyle=:box)

# ╔═╡ 9af435d0-2301-411f-84fc-cc5b151177e1
using Statistics

# ╔═╡ 7ac110de-5ea6-46c9-8a3a-24289b737944
using PlutoUI

# ╔═╡ 763aa97a-bbca-4d89-a5f2-843af83baeb8
using ColorSchemes

# ╔═╡ c61ac63d-969f-4bcb-a0f3-abbafa764923
using AverageShiftedHistograms

# ╔═╡ e7968ebc-8a30-4a9c-b260-fc86883de5a9
using Measures

# ╔═╡ 4d262ca1-e5e9-4009-a445-59395dd25503
using StatsBase

# ╔═╡ 2563b2e6-d2c6-4386-880d-48e8ba31cd44
using LinearAlgebra

# ╔═╡ c57fddaa-5e69-4ea0-8469-1784ab1019bb
using Distributions

# ╔═╡ 0239ffc0-3b91-4441-86ae-cf505742d810
md"""
# Analysis
"""

# ╔═╡ 1b6df2b6-6536-4bec-aa78-39dfb90c7ad6
TableOfContents()

# ╔═╡ a67c7a80-8db6-4d6d-b819-1ee46ebc0396
md"""
# Plotting
"""

# ╔═╡ 5aa567b3-04fa-4aa5-a690-4fb93ccda580
function plot_ore_map(ore_map, cmap=:viridis; kwargs...)
    xl = (1,size(ore_map,1))
    yl = (1,size(ore_map,2))
    return heatmap(ore_map[:,:,1]; title="true ore field", fill=true, clims=(0.0, 1.0), aspect_ratio=1, xlims=xl, ylims=yl, c=cmap, kwargs...)
end

# ╔═╡ 095f3331-54e7-472f-aea9-ab32fe45e680
# function plot_mass_map(m::MineralExplorationPOMDP, ore_map, cmap=:viridis; truth=false, kwargs...)
#     xl = (0.5, size(ore_map,1)+0.5)
#     yl = (0.5, size(ore_map,2)+0.5)
#     s_massive = ore_map .>= m.massive_threshold
#     dim_scale = truth ? 1 : m.dim_scale
#     r_massive = dim_scale*sum(s_massive)
#     mass_fig = heatmap(s_massive[:,:,1], title="massive ore deposits: $(round(r_massive, digits=2))", fill=true, clims=(0.0, 1.0), aspect_ratio=1, xlims=xl, ylims=yl, c=cmap, kwargs...)    
#     return mass_fig
# end

function plot_mass_map(mass_map, r_massive, cmap=:viridis; kwargs...)
    xl = (1,size(mass_map,1))
    yl = (1,size(mass_map,2))
    mass_fig = heatmap(mass_map[:,:,1]; title="massive ore deposits: $(round(r_massive, digits=2))", fill=true, clims=(0.0, 1.0), aspect_ratio=1, xlims=xl, ylims=yl, c=cmap, kwargs...)
    return mass_fig
end

# ╔═╡ 602959ed-946d-48b5-b7d8-ef56e644441c
md"""
# Data
"""

# ╔═╡ f45b80e5-7952-497d-a2af-636ceda95aab
grid_dim = 50

# ╔═╡ df459d63-f23f-48d7-8255-1b3da0910eda
gds = string(grid_dim, "x", grid_dim)

# ╔═╡ 70aad8be-1d90-4d01-bfcd-e98206ce89b0
subdir = "hf_truth_mvb" # "centered_priors_150"

# ╔═╡ 2c82b274-f1ba-4368-b990-03407fe58e02
md"""
## Loading data (blob)
"""

# ╔═╡ 3adf7896-4b2b-4421-aee0-0d19eb64f103
analysis_blob = BSON.load("..\\results\\$subdir\\analysis_blob_$(gds)_pomcp.100.bson")[:analysis];

# ╔═╡ ea904ef3-03fb-4420-b297-c31372624ed4
presults_blob = analysis_blob["results"]

# ╔═╡ d3440665-b3e5-453b-9eb7-41a6859ff4b7
ptiming_blob = map(t->t.time, analysis_blob["timing"])

# ╔═╡ c7a3e631-16bc-4517-bf82-60e15baaa419
mean(first.(presults_blob)), std(first.(presults_blob))

# ╔═╡ d546a273-f8d5-4601-8a1f-4cbbb175ef77
md"""
## Loading data (ellipse)
"""

# ╔═╡ a3b8e58a-660a-433a-9dfb-2c75d8de1635
analysis_ellipse = BSON.load("..\\results\\$subdir\\analysis_ellipse_$(gds)_pomcp.100.bson")[:analysis];

# ╔═╡ 3683782b-a884-44de-8986-66c6d7ba9eba
presults_ellipse = analysis_ellipse["results"]

# ╔═╡ 182672c0-f623-4a4c-80f4-880afaa69b02
ptiming_ellipse = map(t->t.time, analysis_ellipse["timing"])

# ╔═╡ 4296bde1-e654-481c-b365-61d0ddc6dac7
mean(first.(presults_ellipse)), std(first.(presults_ellipse))

# ╔═╡ 3cbac42c-5c87-41e6-99d5-60eeec2530a5
md"""
## Loading data (circle)
"""

# ╔═╡ cd5ea37f-5698-47a8-bb79-82f7c65f455d
analysis_circle = BSON.load("..\\results\\$subdir\\analysis_circle_$(gds)_pomcp.100.bson")[:analysis];

# ╔═╡ c8fd071a-7f74-45c2-aab0-3fd4c0d65dc8
presults_circle = analysis_circle["results"]

# ╔═╡ c3352c24-c8fa-44a4-ab22-3e0edee501bf
ptiming_circle = map(t->t.time, analysis_circle["timing"])

# ╔═╡ 61402b63-a2d7-4d56-8545-aa26ac1a83b3
mean(first.(presults_circle)), std(first.(presults_circle))

# ╔═╡ 3085e225-59b9-4f0c-9586-5d82f307bf34
md"""
## Ore maps
"""

# ╔═╡ 4ec5cb70-6126-479f-b289-7a79fadcaaf5
analysis_selected = analysis_blob

# ╔═╡ 61b8d873-97f9-40d8-8884-78f0b0f03c7b
ore_map = mean(analysis_selected["ore_maps"]);

# ╔═╡ ecd859be-e8db-4380-a6fd-b85438c63781
164.1

# ╔═╡ 8389a09e-bb8b-4946-854b-0f92b58977a8
dim_scale_undo = (0.6^2)

# ╔═╡ f0f3efb8-e52d-4552-9675-762c4f3fedc4
begin
	mass_map = mean(analysis_selected["mass_maps"])
	r_massive = mean(analysis_selected["r_massives"])
end

# ╔═╡ 1192839e-af4e-4ba3-9e3b-a6ebc2fedff6
plot(
	plot_ore_map(ore_map, title="(mean) true ore field"), 
	plot_mass_map(mass_map, r_massive,
		title="(mean) massive ore: $(round(r_massive, digits=2))"),
	size=(700,300))

# ╔═╡ 4e4f6167-875a-41a7-b84a-53d65f5c00a7
md"""
## Relative volume errors
"""

# ╔═╡ 15a60749-135c-4517-af67-3da1cdc62e42
tufte_colors = Dict(
	"pastel_magenta" => colorant"#FF48CF",
	"pastel_purple" => colorant"#8770FE",
	"pastel_blue" => colorant"#1BA1EA",
	"pastel_seagreen" => colorant"#14B57F",
	"pastel_green" => colorant"#3EAA0D",
	"pastel_orange" => colorant"#C38D09",
	"pastel_red" => colorant"#F5615C",	
)

# ╔═╡ 6ca61949-4741-48bb-808a-bf41e6eb1c00
colors = cgrad(:lighttest, 3, categorical=true)
# colors = [
# 	tufte_colors["pastel_blue"],
# 	tufte_colors["pastel_green"],
# 	tufte_colors["pastel_red"],
# ]

# ╔═╡ 7f1e37b7-bdbe-41fe-885b-8abcc0e35578
# :cmyk, :Dark2_3, :darktest, :lighttest

# ╔═╡ c2f7cae6-27ec-4e69-9b3b-410bd63a3443
rel_errs_blob = map(res->res[4], presults_blob);

# ╔═╡ 9a9390c0-4946-40a5-8ed0-a85267ba25f9
rel_errs_ellipse = map(res->res[4], presults_ellipse);

# ╔═╡ 602249bc-4d83-45e3-bd12-cb4097163678
rel_errs_circle = map(res->res[4], presults_circle);

# ╔═╡ 9cbb95db-82ec-4053-9624-3ba34d5cd3e1
Plots.supported_markers()

# ╔═╡ 2e4bc9ef-83c7-45ed-a59e-fec98993e53a
draw_zero() = plot!([xlims()...], [0, 0], lw=1, c=:gray, xlims=xlims(), label=false)

# ╔═╡ 45ad5118-df2e-44b9-9a24-20067944b06c
md"""
## Polar comparison plots
- runtime
- regret
- true positive
- false positive
"""

# ╔═╡ 67db98f5-b5a6-45bd-9676-7b29f8e6f846
r(θ) = 1 + cos(θ) * sin(θ)^2

# ╔═╡ bcfc1720-5737-4145-aec9-cca110544033
plot(r, 0, 2π, lims=(0,1.5), proj=:polar, label=false, xticks=([0, π], ["a","b"]))

# ╔═╡ dab7bfd4-6e8d-4b88-963b-82647f2cd86e
function Plots.gr_polaraxes(rmin::Real, rmax::Real, sp::Plots.Subplot)
    Plots.GR.savestate()
    xaxis = sp[:xaxis]
    yaxis = sp[:yaxis]

    α = 0:45:315
    a = α .+ 90
    sinf = sind.(a)
    cosf = cosd.(a)
    θtick_values, θtick_labels = Plots.get_ticks(sp, xaxis, update = false)
    rtick_values, rtick_labels = Plots.get_ticks(sp, yaxis, update = false)

    #draw angular grid
    if xaxis[:grid]
        Plots.gr_set_line(
            xaxis[:gridlinewidth],
            xaxis[:gridstyle],
            xaxis[:foreground_color_grid],
            sp,
        )
        Plots.gr_set_transparency(xaxis[:foreground_color_grid], xaxis[:gridalpha])
        for i in eachindex(α)
            Plots.GR.polyline([sinf[i], 0], [cosf[i], 0])
        end
    end

    #draw radial grid
    if yaxis[:grid]
        Plots.gr_set_line(
            yaxis[:gridlinewidth],
            yaxis[:gridstyle],
            yaxis[:foreground_color_grid],
            sp,
        )
        Plots.gr_set_transparency(yaxis[:foreground_color_grid], yaxis[:gridalpha])
        for i in eachindex(rtick_values)
            r = (rtick_values[i] - rmin) / (rmax - rmin)
            if r <= 1.0 && r >= 0.0
                Plots.GR.drawarc(-r, r, -r, r, 0, 359)
            end
        end
        Plots.GR.drawarc(-1, 1, -1, 1, 0, 359)
    end

    #prepare to draw ticks
    Plots.gr_set_transparency(1)
    Plots.GR.setlinecolorind(90)
    Plots.GR.settextalign(GR.TEXT_HALIGN_CENTER, GR.TEXT_VALIGN_HALF)

    #draw angular ticks
    if xaxis[:showaxis]
        Plots.GR.drawarc(-1, 1, -1, 1, 0, 359)
        # for i in eachindex(α)
        for i in eachindex(θtick_values)
            θ = θtick_values[i]
            x, y = Plots.GR.wctondc(1.1 * sin(θ), 1.1 * cos(θ))
            # x, y = Plots.GR.wctondc(1.1 * sinf[i], 1.1 * cosf[i])
            # Plots.GR.textext(x, y, string((360 - α[i]) % 360, "^o"))
            Plots.GR.textext(x, y, Plots._cycle(θtick_labels, i))
        end
    end

    #draw radial ticks
    if yaxis[:showaxis]
        for i in eachindex(rtick_values)
            r = (rtick_values[i] - rmin) / (rmax - rmin)
            if r <= 1.0 && r >= 0.0
                x, y = Plots.GR.wctondc(0.05, r)
                Plots.gr_text(x, y, Plots._cycle(rtick_labels, i))
            end
        end
    end
    @info "Moss!"
    Plots.GR.restorestate()
end

# ╔═╡ ce02b5d0-4f40-4dcc-897d-8f49082725af
md"""
## Confusion statistics
"""

# ╔═╡ 499fdc09-fdd6-4ea7-967a-44c4a92a6e4c
function accuracy(presults; extraction_cost=150)
	decisions = last.(presults)
	true_decisions =
		map(res->res[end-1] > extraction_cost ? :mine : :abandon, presults)
	return accuracy(decisions, true_decisions)
end

# ╔═╡ dfa9ae81-5f0a-41d5-8887-712c11a969a4
function count_true_positives(decisions, true_decisions)
	return sum(decisions .== :mine .&& true_decisions .== :mine)
end

# ╔═╡ 43b5f33b-120b-4abe-a970-a21fcf697a7f
function count_true_negatives(decisions, true_decisions)
	return sum(decisions .== :abandon .&& true_decisions .== :abandon)
end

# ╔═╡ 686ebcad-0712-4a2a-ac8f-90389fc5f7b0
function count_false_positives(decisions, true_decisions)
	return sum(decisions .== :mine .&& true_decisions .== :abandon)
end

# ╔═╡ 246dd86d-95bd-42a5-abdb-a7b1280ed42f
function count_false_negatives(decisions, true_decisions)
	return sum(decisions .== :abandon .&& true_decisions .== :mine)
end

# ╔═╡ 6c5bd780-30f1-404c-ac6c-65a9f775636c
function precision(decisions, true_decisions)
	tp = count_true_positives(decisions, true_decisions)
	fp = count_false_positives(decisions, true_decisions)
	return tp / (tp+fp)
end

# ╔═╡ f981ef50-efa6-4716-85e8-ef0a90508d79
function recall(decisions, true_decisions)
	tp = count_true_positives(decisions, true_decisions)
	fn = count_false_negatives(decisions, true_decisions)
	return tp / (tp+fn)
end

# ╔═╡ 9452720c-95d8-43a5-a399-04589b5d5bde
function accuracy(decisions, true_decisions)
	tp = count_true_positives(decisions, true_decisions)
	tn = count_true_negatives(decisions, true_decisions)
	fp = count_false_positives(decisions, true_decisions)
	fn = count_false_negatives(decisions, true_decisions)
	return (tp+tn) / (tp+tn+fp+fn)
end

# ╔═╡ cf93c16a-5547-4815-9adf-ca0147a445e4
md"""
### Confusion matrix
"""

# ╔═╡ 0a4cab3f-131f-4741-80b3-5046390581f6
function confusion_matrix(decisions, true_decisions)
	N = length(decisions)
	tpr = count_true_positives(decisions, true_decisions) / N
	tnr = count_true_negatives(decisions, true_decisions) / N
	fpr = count_false_positives(decisions, true_decisions) / N
	fnr = count_false_negatives(decisions, true_decisions) / N
	return [tpr fpr;
	        fnr tnr]
end

# ╔═╡ dba66cd7-c91d-4619-8ccc-5f39327ed160
function plot_confusion(cm)
	heatmap(rotr90(cm)', c=:viridis, ratio=1, xlims=(0.5, 2.5), ylims=(0.5, 2.5), clims=(0,1))
end

# ╔═╡ 56f9d9a6-82b9-40d4-ba56-fe4768d6b471
md"""
## Plot runtimes
"""

# ╔═╡ 0d57dadb-0681-4fd9-b834-81bfe0042d47
mean(ptiming_blob), std(ptiming_blob)

# ╔═╡ 1bf36d36-e6a8-4aca-8d75-0862d4029acc
mean(ptiming_ellipse), std(ptiming_ellipse)

# ╔═╡ 8951fb50-a5b3-4dce-a5d1-2d79f230f28b
mean(ptiming_circle), std(ptiming_circle)

# ╔═╡ 36476ab4-6c83-49a2-801e-fa442b29f0e1
function plot_runtime(runtime; c=:crimson, label=false, hold=false)
	μ_time = runtime
	σ_time = runtime
	ts = [1:length(μ_time);] .- 1
	plot_fn = hold ? plot! : plot
	plot_fn(ts, μ_time, title="runtime (aggregate)",
		    ribbon=σ_time, fillalpha=0.2,
		    xlabel="rng seed", ylabel="runtime (s)", label=label, lw=2, c=c)
	xlims!(0, ts[end])
	ylims!(minimum(μ_time), ylims()[2])
end

# ╔═╡ fd340989-63ac-4393-ac7b-da9a0d7de25a
begin
	plot_runtime(ptiming_blob)
end

# ╔═╡ 3f4b47b6-77b2-425b-8ce7-fa38f4bde61f
begin
	histogram(ptiming_blob[1:end-1], c=colors[1], alpha=0.5)
	histogram!(ptiming_ellipse[1:end-1], c=colors[2], alpha=0.5)
	histogram!(ptiming_circle[1:end-1], c=colors[3], alpha=0.5)
end

# ╔═╡ 9812856d-eeba-48da-86b4-b8005a577d5e
begin
	m_ash=40
	μσ_label = d -> string(round(mean(d), digits=2), " ± ", round(std(d), digits=2))
	plot_fit = (d,c,lab,hold=false) -> begin
		plot_fn = hold ? plot! : plot
		plot_fn(ash(d, m=m_ash), hist=false, label="$lab ($(μσ_label(d)))", c=c)
	end
	plot_fit(ptiming_blob, colors[1], "blob")
	plot_fit(ptiming_ellipse, colors[2], "ellipse", true)
	plot_fit(ptiming_circle, colors[3], "circle", true)
end

# ╔═╡ 94f09ea2-792a-4936-8f40-292d4fbd376b
begin
	bar(["blob" "ellipse" "circle"],
		[mean(ptiming_blob) mean(ptiming_ellipse) mean(ptiming_circle)],
		color=reshape(collect(colors), (1, 3)))
end

# ╔═╡ e9854a55-31ad-4d7f-8331-e5575a6375c6
md"""
# 3D value plot
"""

# ╔═╡ 8333262e-8968-4396-8dbf-7a30c5d2f466
md"""
## Load sweep data (blob/blob)
"""

# ╔═╡ 1b514edc-cde4-4b1f-add8-b21ef545b494
# analysis_sweep_blob_blob = BSON.load("..\\src\\analysis_blob_blob_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

# ╔═╡ 5002c00d-c831-4d9a-8c4a-5aaa16ac8819
# plot_value_grid(analysis_sweep_blob_blob, "blob/blob")

# ╔═╡ f588a53c-3b72-4f98-89ac-219555d6791c
md"""
## Load sweep data (blob/ellipse)
"""

# ╔═╡ eb14649a-b273-469c-a720-02480641283a
# analysis_sweep_blob_ellipse = BSON.load("..\\src\\analysis_blob_ellipse_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

# ╔═╡ 7acb3d44-0a71-4a5f-9678-efc716f530d0
# plot_value_grid(analysis_sweep_blob_ellipse, "blob/ellipse")

# ╔═╡ f85ead75-cc42-4ed5-a410-7c51f6e8a6a5
md"""
## Load sweep data (blob/circle)
"""

# ╔═╡ 462df3d0-860f-4629-af3e-95648f56c1df
# analysis_sweep_blob_circle = BSON.load("..\\src\\analysis_blob_circle_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

# ╔═╡ b10e3cb2-f8f1-4743-984c-4dc167aa7f59
# plot_value_grid(analysis_sweep_blob_circle, "blob/circle")

# ╔═╡ e2052af4-26d8-45ed-aeed-0ab0dd37b6f4
md"""
## Sweep matrix (truth and belief)
"""

# ╔═╡ d9b4cbe2-884c-407f-843c-f5cbda9c3b14
combined_results = [(1,"one","1"), (2, "two", "2"), (3, "three", "3")]

# ╔═╡ 97ff0196-af8e-4060-a590-2784e3a967b1
[[combined_results[i][j] for i in 1:length(combined_results)] for j in 1:length(combined_results[1])]

# ╔═╡ 4f04be51-8452-4584-b1c8-804e957de851
# begin
# 	analysis_sweep_ellipse_blob = BSON.load("..\\src\\analysis_ellipse_blob_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

# 	analysis_sweep_ellipse_ellipse = BSON.load("..\\src\\analysis_ellipse_ellipse_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

# 	analysis_sweep_ellipse_circle = BSON.load("..\\src\\analysis_ellipse_circle_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];


# 	analysis_sweep_circle_blob = BSON.load("..\\src\\analysis_circle_blob_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

# 	analysis_sweep_circle_ellipse = BSON.load("..\\src\\analysis_circle_ellipse_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

# 	analysis_sweep_circle_circle = BSON.load("..\\src\\analysis_circle_circle_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];
# end;

# ╔═╡ 6570b47d-e84a-4c63-a0d1-a831b647d8d3
# begin
# 	mfunc = (data,ttl) -> begin
# 		plot_value_grid(data, ttl;
# 			reduced=true, cmap=shared_cmap, clims=(svmin,svmax))
# 	end

#     mplot = plot(
#         mfunc(analysis_sweep_blob_blob, "blob/blob"),
#         mfunc(analysis_sweep_blob_ellipse, "blob/ellipse"),
#         mfunc(analysis_sweep_blob_circle, "blob/circle"),
    
#         mfunc(analysis_sweep_ellipse_blob, "ellipse/blob"),
#         mfunc(analysis_sweep_ellipse_ellipse, "ellipse/ellipse"),
#         mfunc(analysis_sweep_ellipse_circle, "ellipse/circle"),
    
#         mfunc(analysis_sweep_circle_blob, "circle/blob"),
#         mfunc(analysis_sweep_circle_ellipse, "circle/ellipse"),
#         mfunc(analysis_sweep_circle_circle, "circle/circle"),
    
#         layout=@layout[A B C; D E F; G H I], size=(700,700)
#     )
# end

# ╔═╡ ce739257-3473-426f-aef4-86a7e016a575
md"""
## Bayesian analysis using beta distributions
"""

# ╔═╡ cfe230ef-cf77-4a45-9c6b-85781383203d
begin
    data5_blob_blob = BSON.load("..\\results\\$subdir\\analysis_blob_blob_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

    data5_blob_ellipse = BSON.load("..\\results\\$subdir\\analysis_blob_ellipse_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

    data5_blob_circle = BSON.load("..\\results\\$subdir\\analysis_blob_circle_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];


    # data5_ellipse_blob = BSON.load("..\\src\\analysis_ellipse_blob_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

    # data5_ellipse_ellipse = BSON.load("..\\src\\analysis_ellipse_ellipse_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

    # data5_ellipse_circle = BSON.load("..\\src\\analysis_ellipse_circle_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];


    # data5_circle_blob = BSON.load("..\\src\\analysis_circle_blob_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

    # data5_circle_ellipse = BSON.load("..\\src\\analysis_circle_ellipse_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];

    # data5_circle_circle = BSON.load("..\\src\\analysis_circle_circle_seeds5_sweep_xy.10.30.50_pomcp.10.100.1000.bson")[:analysis];
end;

# ╔═╡ 232b0feb-c67d-4ccf-ad2b-bba1802e9344
data5 = [
	data5_blob_blob,
	data5_blob_ellipse,
	data5_blob_circle,
	# data5_ellipse_blob,
	# data5_ellipse_ellipse,
	# data5_ellipse_circle,
	# data5_circle_blob,
	# data5_circle_ellipse,
	# data5_circle_circle,
]

# ╔═╡ d972f53b-90fc-4d12-8697-1b18db694a62
data5[1]

# ╔═╡ cf562730-8ac6-4b45-a311-a8208c3982fb
titles = [
    "blob/blob",
    "blob/ellipse",
    "blob/circle",
    # "ellipse/blob",
    # "ellipse/ellipse",
    # "ellipse/circle",
    # "circle/blob",
    # "circle/ellipse",
    # "circle/circle",
]

# ╔═╡ 78d81469-4c41-4be9-afbc-c5ce52e9b2db
function uniform_betas(data5)
	return [[Beta(1,1) for i in 1:length(d["grid_dim_xys"])] for d in data5]
end

# ╔═╡ cf75c166-f776-4b49-b3d2-840aaf8b2231
betas = uniform_betas(data5)

# ╔═╡ fb8ba111-ce93-4c73-8602-4b6199eba3e0
function findindexmulti(data, grid_dim_xy, pomcpow_iter)
	for i in 1:length(data["grid_dim_xys"])
		if data["grid_dim_xys"][i][1] == grid_dim_xy &&
				data["pomcpow_iters"][i][1] == pomcpow_iter
			return i
		end
	end
	return missing
end

# ╔═╡ 598875f9-e08f-440c-907b-8b840080988c
begin
	data = data5[1]
	extraction_cost = 150 # TODO: from m.extraction_cost
	massive_and_decision = 
		map(output->output[end-1:end], data["results"][findindexmulti(data, 50, 1000)])
	decisions = last.(massive_and_decision)
	true_decisions =
		map(md->md[1] > extraction_cost ? :mine : :abandon, massive_and_decision)
end

# ╔═╡ f5b47f3e-fa1e-4c33-9947-14b8365e293e
count_true_positives(decisions, true_decisions)

# ╔═╡ 29de6acb-cb1a-42d7-ae8d-cfeb09be2136
count_true_negatives(decisions, true_decisions)

# ╔═╡ f68f7701-746f-4ffa-9a80-6b6d12f85559
count_false_positives(decisions, true_decisions)

# ╔═╡ 869bd9fa-3c23-4d0a-b7fa-7e5942b0a556
count_false_negatives(decisions, true_decisions)

# ╔═╡ 0a458ee4-6781-4f79-94c4-3d5d0d698e9b
precision(decisions, true_decisions)

# ╔═╡ 6c495b97-74b6-4dd5-a66d-2b626899ea2b
recall(decisions, true_decisions)

# ╔═╡ 3318952d-9a71-4f7c-bf84-7e3bbf9e1983
massive_and_decision

# ╔═╡ 7229656e-3525-4be4-a11f-dd277c3b46cd
cm = confusion_matrix(decisions, true_decisions)

# ╔═╡ 9a02c71a-b8cf-41fc-bca7-490b2c2f8966
function accuracy(data, x, y, s=missing; extraction_cost=150)
	i = findindexmulti(data, x, y)
	massive_and_decision = 
		map(output->output[end-1:end], data["results"][i])
	decisions = last.(massive_and_decision)
	true_decisions =
		map(md->md[1] > extraction_cost ? :mine : :abandon, massive_and_decision)
	return accuracy(decisions, true_decisions)
end

# ╔═╡ f0411d8c-a907-4ca7-9d49-a965c884aaff
accuracy(decisions, true_decisions)

# ╔═╡ 1b9633f7-f2f6-4de8-9af9-cc325da589fb
accuracy(data5[1], 10, 10)

# ╔═╡ 53e84c72-293a-4818-853a-a9da2e7e580d
function betavalue(data, x, y, idx_shapes; true_decision=:mine)
	i = findindexmulti(data, x, y)
	s = idx_shapes
	return mean(betas[s][i]) # Bayesian
	# return (betas[s][i].α - 1) / (betas[s][i].α + betas[s][i].β - 2) # Frequentist
end

# ╔═╡ cc0afc27-926d-4fc7-b802-ea31faaa7041
betavalue(data5[3], 30, 100, 3)

# ╔═╡ 7a7bc902-5439-4248-9fc3-3baaacab4531
function update_beta!(data, x, y, idx_shapes; true_decision=:mine)
	global betas
	i = findindexmulti(data, x, y)
	decisions = map(d->d[end], data["results"][i])
	s = idx_shapes
	for decision in decisions
		r = decision == true_decision
		betas[s][i] = Beta(betas[s][i].α + r, betas[s][i].β + 1 - r)
	end
	return betas
end

# ╔═╡ 272b238e-c9f6-4856-8d69-1ec8d07c1b0b
begin
	for x in [10, 30, 50]
		for y in [10, 100, 1000]
			map(s->update_beta!(data5[s], x, y, s), 1:length(data5))
		end
	end
	betas[3]
end

# ╔═╡ 59b6f645-616a-4694-9c43-915926d9bf74
function get_data_xy(data)
	X = unique(data["grid_dim_xys"])
	Y = unique(data["pomcpow_iters"])

	X = mapreduce(unique, vcat, X)
	Y = mapreduce(unique, vcat, Y)
	return X, Y
end

# ╔═╡ 83a79501-d3ba-4adf-b183-4d1b28b815bd
md"""
## Plotting beta means
"""

# ╔═╡ fcf97cb8-b9c4-40f6-a694-af848046d4bd


# ╔═╡ 1e3318e5-3181-4d8a-8602-22291b3750ef
# begin
# 	mbfunc = s -> begin
# 		plot_betas(data5[s], titles[s], s; reduced=true)
# 	end

# 	mbplot = plot([mbfunc(s) for s in 1:length(data5)]...,    
#         layout=@layout[A B C; D E F; G H I], size=(700,700)
#     )
# 	pcbar = scatter([0], [0], alpha=0, 
# 		zcolor=[0,1], clims=(0,1), c=get_colorbar(0, 1; vmid=0.5),
# 		axis=false, tick=nothing, label=false)
# 	pylabel = scatter([0], [0], alpha=0,
# 		axis=false, tick=nothing, label=false,
# 	    ylabel="planning iterations")
# 	pxlabel = scatter([0], [0], alpha=0,
# 		axis=false, tick=nothing, label=false,
# 	    xlabel="grid dimension (n×n)")
# 	plot(pylabel, mbplot, pcbar, pxlabel,
# 		 layout=@layout([a{0.00001w} b c{0.1w}; d{0.01h}]),
#  		 size=(700,700), plot_title="expected true positive rate [truth/belief]")
# end

# ╔═╡ 53d2a31e-0fd5-49a9-ae81-5a55b22e328d
md"""
## Regret

$$\begin{gather}
R_\text{best} = \max\left\{0, \text{massive ore} - \text{extraction cost}\right\}\\
\operatorname{regret} = R_\text{best} - R
\end{gather}$$
"""

# ╔═╡ ac753e4d-e7f6-4eb6-af22-d06f58473fb7
function regret(data, x, y)
	i = findindexmulti(data, x, y)
	results = data["results"]
	extraction_cost = 150.0
	massive_ore = map(res->res[end-1], results[i])
	r_best = max.(0, massive_ore .- extraction_cost)
	returns = map(res->res[1], results[i])
	return r_best .- returns
end

# ╔═╡ d1188bf0-92c1-4790-a6a6-953a18f46354
regrets = regret(data5[3], 10, 100)

# ╔═╡ abd7e135-f4b5-4cbe-86da-19e94039a7e7
mean(regrets)

# ╔═╡ c623ed18-c95c-4ef3-88c8-57a1213660bf
quantile(regrets)

# ╔═╡ 4a552377-db94-47df-a438-2309bcb7b814
data5[1]["results"][findindexmulti(data5[1], 10, 10)]

# ╔═╡ 74fdc433-91ff-4c08-9c46-2563be98f317
regret_fn_90th_quant = (d,x,y)->quantile(regret(d,x,y), 0.9)

# ╔═╡ c1bccfa7-9f45-4ed2-8a4b-5a269329ed0a
quantile(regret(data5[3], 50, 10), 0.2)

# ╔═╡ 04144789-7339-4f3e-943b-957d87e6dc33
# function get_colorbar_regret(data::Dict, X::Vector, Y::Vector)
function get_colorbar_regret(x, y)
	# Z = [betavalue(data, x, y, s) for x in X for y in Y]
    # vmin = minimum(Z)
    # vmax = maximum(Z)
    # buckets = [vmin, vmin/2, vmid, vmax/2, vmax] # shift colormap so 0 is at center
    # normed = (buckets .- vmin) / (vmax - vmin)
	# return cgrad(:jet, normed, rev=true)
	# return cgrad([:darkred, :red, :white, :forestgreen, :black], normed)
	return cgrad(:jet)
end

# ╔═╡ 7017212c-71da-4fe2-90c6-fe95534446da
function plot_regret(data, title, idx_shapes;
		reduced=false, cmap=nothing, clims=nothing,
		show_cbar=!reduced, show_xlabel=!reduced, show_ylabel=!reduced,
		value_func=(d,x,y)->mean(regret(d,x,y)))
	plot()
	if reduced
		title!(title)
	else
		title!("expected true positive rate [$title]")
	end

	if show_xlabel
		xlabel!("grid dimension (n×n)")
	end
	if show_ylabel
		ylabel!("planning iterations")
	end
	
	X, Y = get_data_xy(data)
	
	# Color scheme (centered at 0)
	if isnothing(cmap)
		# Based off individual values of the data (i.e., not shared)
		cmap = get_colorbar_regret(nothing, nothing)
	end
	
	contourf!(X, Y, (x,y)->value_func(data, x, y),
		      c=cmap, cbar=show_cbar, clims=clims, levels=16,
			  size=(500,400), yaxis=:log, linewidth=0.25, linecolor=:black)

	saved_lims = (xlims(), ylims())
	G = [(x, y) for x in X, y in Y]
	_X, _Y = first.(G), last.(G)
	scatter!(_X, _Y, ms=4, color=:white, label=false, marker=:circle)
	
	if reduced
		xlims!(saved_lims[1]...)
		ylims!(saved_lims[2]...)
	end

	return plot!()
end

# ╔═╡ c21fe459-c5c6-4324-ac17-005dddc94bb7
md"""
# Pareto curves (regret vs. runtime)
"""

# ╔═╡ d06f41df-45ee-4ed0-adb1-ddd5666bd243
[map(t->t.time, timing) for timing in data5[1]["timing"]]

# ╔═╡ d79d7b50-6714-4175-8eba-ef3ad7451de9
data5[1]

# ╔═╡ c69f4e20-0489-4d93-98e0-e47a78a9fb55
function goal_programming(𝐱::Dict, 𝐲::Dict, goal=[0,0]; p=2)
	return argmin(kv->norm(last(kv) - goal, p), merge(vcat, 𝐱, 𝐲))
end

# ╔═╡ e3d813b2-1c52-4d05-8953-50ed8412cc6e
begin
	plot()
	mean_timings = Dict()
	std_timings = Dict()
	mean_regrets = Dict()
	std_regrets = Dict()
	is_mins = true
	std_scale = 10
	for s in [1,2,3]
		for xy in [10, 30, 50]
			for iters in [10, 100, 1000]
				i = findindexmulti(data5[s], xy, iters)
				time_divisor = is_mins ? 60 : 1
				times = map(t->t.time/time_divisor, data5[s]["timing"][i])
				mean_time = mean(times)
				std_time = std(times)

				regrets = regret(data5[s], xy, iters)
				mean_regret = mean(regrets)
				std_regret = std(regrets)
				
				key = "$s, $(xy)x$(xy), $iters"
				mean_timings[key] = mean_time
				mean_regrets[key] = mean_regret
				std_timings[key] = std_time
				std_regrets[key] = std_regret
				# @info "$(xy)x$(xy), $iters, $mean_time"
				if s == 1
					mark = :star6
				elseif s == 2
					mark = :diamond
				elseif s == 3
					mark = :circle
				end
				if xy == 10
					c = :red
				elseif xy == 30
					c = :green
				elseif xy == 50
					c = :blue
				end
				if iters == 10
					ms = 4
				elseif iters == 100
					ms = 8
				elseif iters == 1000
					ms = 16
				end
				label = iters == 100 && s == 3 ? "$(xy)x$(xy)" : false
				scatter!([mean_time], [mean_regret],
						 xerr=std_time/std_scale, yerr=std_regret/std_scale, msc=c,
					     c=c, ms=ms, alpha=0.3, label=label, marker=mark)
			end
		end
	end
	# mean_timings, mean_regrets
	pareto_optimal = last(goal_programming(mean_timings, mean_regrets; p=Inf))
	plot!([0, pareto_optimal[1]], [0, pareto_optimal[2]],
		   c=:black, lw=1, style=:dash, label="pareto optimal")
	plot!(legend=:topright)
	title!("pareto curve (error σ/$std_scale)")
	xlabel!("mean runtime ($(is_mins ? "mins." : "secs."))")
	ylabel!("mean regret")
	xlims!(0,xlims()[end])
	ylims!(0,ylims()[end])
	# hline!([0], c=:gray, label=false)
	# vline!([0], c=:gray, label=false)
end

# ╔═╡ bf75df30-32cf-49c0-8048-b891f5f03e13
last(pareto_optimal)

# ╔═╡ 94a167e5-7a24-4acc-8fdf-3fc5ad3f4af7
# begin
# 	plot()
# 	for xy in [10, 30, 50]
# 		for iters in [10, 100, 1000]
# 			key = "$(xy)x$(xy), $iters"
# 			mean_time = mean_timings[key]
# 			mean_regret = mean_regrets[key]
# 			if xy == 10
# 				c = :red
# 			elseif xy == 30
# 				c = :green
# 			elseif xy == 50
# 				c = :blue
# 			end
# 			if iters == 10
# 				ms = 2
# 			elseif iters == 100
# 				ms = 4
# 			elseif iters == 1000
# 				ms = 8
# 			end
# 			label = iters == 100 ? "$(xy)x$(xy)" : false
# 			scatter!([mean_time], [mean_regret], c=c, ms=ms, alpha=0.5, label=label)
# 		end
# 	end
# 	plot!(legend=:bottomright)
# 	title!("pareto curve [blob/circle]")
# 	xlabel!("mean regret")
# 	ylabel!("mean runtime")
# 	xlims!(-10,xlims()[end])
# 	ylims!(-10,ylims()[end])
# 	hline!([0], c=:gray, label=false)
# 	vline!([0], c=:gray, label=false)
# end

# ╔═╡ b7d46821-abb0-42ba-88a1-f9e24371a455
md"""
# Reward function testing
"""

# ╔═╡ 207d023f-bd15-4e73-9de8-6291dea3ec39
@bind abs_err Slider(0:0.1:30, default=0, show_value=true)

# ╔═╡ 39803ea6-ecd5-48c4-a4b7-729ab2743263
function reward2(answer, abs_err)
	# x = sign(answer) * abs_err
	# return x / (1 + exp(-x))
	return -abs_err^2
end

# ╔═╡ 600eae30-53bf-46b1-863e-ca28b64d1834
reward2(+1, abs_err)

# ╔═╡ 4f14d88e-7ef1-4ab2-992b-5afbaddb3e37
reward2(-1, abs_err)

# ╔═╡ 367369f3-c2cb-4c56-b4fe-e29c9536289f
begin
	reward_X = -10:0.1:10
	plot(reward_X, x->reward2(sign(x), abs(x)),
		 c=:black, lw=2, xlabel="abs. error", ylabel="reward", label=false,
	)
	draw_zero()
end

# ╔═╡ 0b91ab7d-f5bf-4889-a7ee-d1b5783e72e2
_beta = Beta(1,1)

# ╔═╡ 914a81ee-e8c5-4173-bceb-710dbfedc4a3
mean(_beta)

# ╔═╡ 08102844-7ba0-4ca2-97fc-7b66410859f4
plot(0:0.01:1, x->pdf(Beta(1,1), x))

# ╔═╡ d8723841-f4b2-4c09-bac3-3295b8e5176d
function reward(answer, abs_err)
	# return sign(answer)
	# if sign(answer) == 1
	# 	return mean(Beta(2,1))
	# else
	# 	return mean(Beta(1,2))
	# end
	return sign(answer) * abs_err
	# return sign(answer) * abs_err^2
	# if abs_err == 0
	# 	10
	# elseif sign(answer) == 1
	# 	100inv(abs_err)
	# else
	# 	-sqrt(abs_err)
	# end
	# sign(answer) * (abs_err == 0 ? 10 : inv(abs_err))
end

# ╔═╡ 4d1911a1-2c79-4524-a630-5eaf9e462ffe
begin
	plot(reward_X, x->reward(sign(x), abs(x)),
		 c=:black, lw=2, xlabel="abs. error", ylabel="reward", label=false,
	)
	draw_zero()
end

# ╔═╡ 166ad666-5a03-410c-9448-30c518051d3b
md"""
$$R(s, a) = \begin{cases}
x & \text{if}\\
y & \text{otherwise}
\end{cases}$$
"""

# ╔═╡ 7b56bf68-6d3c-4b6e-8826-082215409208
md"""
## Sweep utilities
"""

# ╔═╡ 172f1589-5a7c-4dd9-aac4-3be044dd32df
datasets = [
	analysis_sweep_blob_blob,
	analysis_sweep_blob_ellipse,
	analysis_sweep_blob_circle,
	analysis_sweep_ellipse_blob,
	analysis_sweep_ellipse_ellipse,
	analysis_sweep_ellipse_circle,
	analysis_sweep_circle_blob,
	analysis_sweep_circle_ellipse,
	analysis_sweep_circle_circle,
];

# ╔═╡ 73abd22c-7438-4398-bc43-6441ca83939f
function get_colorbar(vmin::Real, vmax::Real; vmid=0)
    buckets = [vmin, vmin/2, vmid, vmax/2, vmax] # shift colormap so 0 is at center
    normed = (buckets .- vmin) / (vmax - vmin)
    # cmap = cgrad(:RdYlGn, normed)
	# return cgrad([:darkred, :red, :white, :green, :darkgreen], normed)

	# return cgrad(:jet, normed, rev=true)
	return cgrad([:darkred, :red, :white, :forestgreen, :black], normed)
end

# ╔═╡ 45da8856-59b3-4dab-9d9f-edae3269baad
function get_colorbar(data::Dict)
	X = unique(data["grid_dim_xys"])
	Y = unique(data["pomcpow_iters"])
	return get_colorbar(data, X, Y)
end

# ╔═╡ 3afb0730-51cb-4f7f-a26b-8618df4126d8
data5[1]

# ╔═╡ 19236596-7225-4767-8348-4be4caa684c2
mfunc(analysis_sweep_ellipse_circle, "ellipse/circle")

# ╔═╡ 96182d6b-4821-4abc-8a41-3be25b9f6e0b
function findindex(data, grid_dim_xy, pomcpow_iter)
	for i in 1:length(data["grid_dim_xys"])
		if data["grid_dim_xys"][i] == grid_dim_xy &&
				data["pomcpow_iters"][i] == pomcpow_iter
			return i
		end
	end
	return missing
end

# ╔═╡ 0e39f099-62d3-49ba-814d-9df55be91f81
function value(data, x, y; true_decision=:mine)
	i = findindex(data, x, y)
	rel_err = data["results"][i][4][end]
	abs_err = data["results"][i][3][end]
	# return rel_err
	decision = data["results"][i][end]
	answer = decision == true_decision ? 1 : -1
	# return answer ? abs_err : -abs_err
	return reward(answer, abs_err)
	# return rel_err
	# return decision == true_decision ? 1 : -1
end

# ╔═╡ 70b883ec-15b7-432c-b375-d5519ea32028
function get_colorbar(data::Dict, X::Vector, Y::Vector)
	Z = [value(data, x, y) for x in X for y in Y]
    vmin = minimum(Z)
    vmax = maximum(Z)
	return get_colorbar(vmin, vmax)
end

# ╔═╡ 5c84a2ac-b719-47fd-b406-81a2255aa9c2
function get_colorbar_beta(data::Dict, X::Vector, Y::Vector, s)
	# Z = [betavalue(data, x, y, s) for x in X for y in Y]
    # vmin = minimum(Z)
    # vmax = maximum(Z)
	return get_colorbar(0, 1; vmid=0.5)
end

# ╔═╡ a5c975e4-8118-421c-a55e-a4d438337c71
function plot_confusion(data, subtitle, idx_shapes;
		reduced=false, cmap=nothing,
		show_cbar=!reduced, show_xlabel=!reduced, show_ylabel=!reduced,
		value_fn=accuracy,
		title="accuracy",
		# betavalue,
		# truepositivevalue, truenegativevalue, falsepositivevalue, falsenegativevalue
	)
	plot()
	if reduced
		title!(subtitle)
	else
		title!("$title [$subtitle]")
	end

	if show_xlabel
		xlabel!("grid dimension (n×n)")
	end
	if show_ylabel
		ylabel!("planning iterations")
	end
	
	X, Y = get_data_xy(data)
	
	# Color scheme (centered at 0)
	if isnothing(cmap)
		# Based off individual values of the data (i.e., not shared)
		cmap = get_colorbar_beta(data, X, Y, idx_shapes)
	end
	
	contourf!(X, Y, (x,y)->value_fn(data, x, y, idx_shapes),
		      c=cmap, cbar=show_cbar, clims=(0,1), levels=16,
			  size=(500,400), yaxis=:log, linewidth=0.25, linecolor=:black)

	saved_lims = (xlims(), ylims())
	G = [(x, y) for x in X, y in Y]
	_X, _Y = first.(G), last.(G)
	scatter!(_X, _Y, ms=4, color=:white, label=false, marker=:circle)
	
	if reduced
		xlims!(saved_lims[1]...)
		ylims!(saved_lims[2]...)
	end

	return plot!()
end

# ╔═╡ 02c3f953-69b7-4419-96d2-18341d8f47f8
begin
	plot_confusion(cm)
	title!("10x10 | 100 iters")
end

# ╔═╡ a3ba5f18-7d8e-4a57-85e1-56933771bd77
begin
	s = 1
	plot_confusion(data5[s], titles[s], s)
end

# ╔═╡ 94364d2e-7936-4919-9641-5d1ec1762f36
data5[s]

# ╔═╡ dc4e45cc-eb23-4dc5-83dc-a2ab136db691
[begin
	X, Y = get_data_xy(data5[s])
	G = [(x, y) for y in Y, x in X]
	(mean(b), G[i])
end for (i,b) in enumerate(betas[s])]

# ╔═╡ 989dea24-8610-4275-b4fb-a601a4bc870d
begin
	mbfunc = (s;kwargs...) -> begin
		plot_confusion(data5[s], titles[s], s; reduced=true, kwargs...)
	end

	mbtitle = plot(title="accuracy [truth/belief]",
		           grid=false, axis=false, tick=nothing, bottom_margin=-50Plots.px)

	mbplot = plot([mbfunc(s; show_ylabel=(s==1), show_xlabel=(s==2))
		           for s in 1:length(data5)]...,
        layout=@layout[A B C], size=(700,210)
    )
	pcbar = scatter([0], [0], alpha=0, 
		zcolor=[0,1], clims=(0,1), c=get_colorbar(0, 1; vmid=0.5),
		axis=false, tick=nothing, label=false)
	plot(mbtitle, mbplot, pcbar,
		 layout=@layout([a{0.01h}; b c{0.1w}]),
 		 size=(710,280), bottom_margin=5mm)
end

# ╔═╡ 9d7caa02-54a1-40f0-816a-85faff5992d6
function bounded_values(datasets; value_fn=(d,x,y)->value(d,x,y))
	bounds = []
	for m_fn in [minimum, maximum]
		bound = m_fn(data->begin
			X = unique(data["grid_dim_xys"])
			Y = unique(data["pomcpow_iters"])
			m_fn(value_fn(data, x, y) for x in X for y in Y)
		end, datasets)
		push!(bounds, bound)
	end
	return bounds
end

# ╔═╡ 0ad9c8b0-918f-4e4c-b9a9-b41e296d7dd4
function get_colorbar_bounded(datasets::Vector; value_fn)
	vmin, vmax = bounded_values(datasets; value_fn=value_fn)
	return get_colorbar(vmin, vmax), vmin, vmax
end

# ╔═╡ ba97c361-013a-41f3-a54c-02460ef63270
shared_cmap, svmin, svmax = get_colorbar_bounded(datasets; value_fn=value);shared_cmap

# ╔═╡ 10e845f2-2551-4453-a7c2-180118503283
svmin, svmax

# ╔═╡ 532beb9a-06c0-4ef7-b905-813aaec5a55d
function bounded_values_multi(datasets; value_fn=(d,x,y)->value(d,x,y))
	bounds = []
	for m_fn in [minimum, maximum]
		bound = m_fn(data->begin
			X, Y = get_data_xy(data)
			m_fn(value_fn(data, x, y) for x in X for y in Y)
		end, datasets)
		push!(bounds, bound)
	end
	return bounds
end

# ╔═╡ ef812da1-135c-4664-b79b-91430f46d3dc
function get_colorbar_bounded_multi(datasets::Vector; value_fn, lower=nothing)
	vmin, vmax = bounded_values_multi(datasets; value_fn=value_fn)
	if !isnothing(lower)
		vmin = lower # generally for cases when we want a shared 0 lower bound
	end
	return get_colorbar_regret(vmin, vmax), vmin, vmax
end

# ╔═╡ ea29156a-4e96-4f7c-a9ca-af235259d14b
function plot_sweep_regret(dataset, regret_fn, regret_title; colorbar_regret_fn=regret_fn, kwargs...)
	shared_cmap_regret, svmin_reget, svmax_regret =	get_colorbar_bounded_multi(data5; value_fn=colorbar_regret_fn, lower=0)

	regret_clims = (svmin_reget,svmax_regret)
	rgfunc = (s;kwargs...) -> begin
		plot_regret(dataset[s], titles[s], s; reduced=true, value_func=regret_fn,
			cmap=shared_cmap_regret, clims=regret_clims, kwargs...)
	end

	rgtitle = plot(title="$regret_title [truth/belief]",
		           grid=false, axis=false, tick=nothing, bottom_margin=-50Plots.px)

	rgplot = plot([rgfunc(s; show_ylabel=(s==1), show_xlabel=(s==2))
		           for s in 1:length(dataset)]...,
        layout=@layout[A B C])
	rgcbar = scatter([0], [0], alpha=0, 
		zcolor=[regret_clims...], clims=regret_clims, c=shared_cmap_regret,
		axis=false, tick=nothing, label=false)
	plot(rgtitle, rgplot, rgcbar,
		 layout=@layout([a{0.01h}; b c{0.1w}]),
 		 size=(710,280), bottom_margin=5mm)
end

# ╔═╡ b3296ce7-f606-4c43-b10a-3603c887be4e
begin
	regret_fn_mean = (d,x,y)->mean(regret(d,x,y))
	regret_title_mean = "expected regret"
	p_regret = plot_sweep_regret(data5, regret_fn_mean, regret_title_mean;
		colorbar_regret_fn=regret_fn_90th_quant)
end

# ╔═╡ eaebf4df-c915-4136-aa99-7101563d4aca
plot_sweep_regret(data5, regret_fn_mean, regret_title_mean)

# ╔═╡ 94b106f6-409d-457f-a4a7-4ad68a9fd0ef
begin
	α_quantile = 0.8 # 0.8? 0.2?
	regret_fn_quant = (d,x,y)->quantile(regret(d,x,y), α_quantile)
	regret_title_quant = "$(round(Int, 100*α_quantile))th percentile of regret"
	p_regret80 = plot_sweep_regret(data5, regret_fn_quant, regret_title_quant;
		colorbar_regret_fn=regret_fn_90th_quant)
end

# ╔═╡ d76d272b-d810-4726-86ab-fa89f81bcb12
quantile(regrets, α_quantile)

# ╔═╡ f57dde29-d588-4fa5-bd41-f9b13b9d877e
begin
	plot(x->ecdf(regrets)(x), label=false, fill=true, 0, 150, size=(300,200), fillcolor=:gray, c=:black, title="cumulative regret (50×50, 10 iters.)", titlefontsize =10)
	cr_xlims = xlims()
	cr_ylims = ylims()
	low_quant_regret = quantile(regrets, α_quantile)
	plot!([cr_xlims[1], low_quant_regret], [α_quantile, α_quantile], lw=4, c=:crimson, lab=false)
	plot!([low_quant_regret,low_quant_regret], [0, α_quantile], lw=4, c=:crimson, lab=false)
	scatter!([low_quant_regret], [α_quantile], ms=4, c=:red, lab=false)
	xlims!(cr_xlims...)
	ylims!(0, cr_ylims[2])
end

# ╔═╡ 3bad29de-edc2-4248-b5fb-5c1aa360c457
begin
	α_quantile2 = 0.9 # 0.8? 0.2?
	regret_fn_quant90 = (d,x,y)->quantile(regret(d,x,y), α_quantile2)
	regret_title_quant90 = "$(round(Int, 100*α_quantile2))th percentile of regret"
	p_regret90 = plot_sweep_regret(data5, regret_fn_quant90, regret_title_quant90;
		colorbar_regret_fn=regret_fn_90th_quant)
end

# ╔═╡ 9e169ad0-b323-4cd8-a60f-fbbd35bca36a
bounded_values_multi(data5; value_fn=(d,x,y)->mean(regret(d,x,y)))

# ╔═╡ 5ad33171-5466-4165-9615-6b618eb4f5f8
function plot_value_grid(data, title; reduced=false, cmap=nothing, clims=nothing)
	plot()
	if reduced
		title!(title)
	else
		title!("sign(true positive) × abs. error [$title]")
		xlabel!("grid dimension (n×n)")
		ylabel!("planning iterations")
	end
	# p_curr, v_curr = s

	X = unique(data["grid_dim_xys"])
	Y = unique(data["pomcpow_iters"])
	
	# Color scheme (centered at 0)
	if isnothing(cmap)
		# Based off individual values of the data (i.e., not shared)
		cmap = get_colorbar(data, X, Y)
	end
	
	# grid.cutPoints...
	contourf!(X, Y, (x,y)->value(data, x,y),
		      c=cmap, cbar=!reduced, clims=clims,
			  size=(500,400), yaxis=:log, linewidth=0.25, linecolor=:black)

	# idx_curr, weights_curr = interpolants(grid, s)
	# for (i,v) in enumerate(idx_curr)
	# 	xₚ, yₚ = grid[v]
	# 	w = weights_curr[i]
	# 	plot!([p_curr, xₚ], [v_curr, yₚ], lw=10w, label=false, color=get(gradient, w))
	# end

	saved_lims = (xlims(), ylims())
	G = [(x, y) for x in X, y in Y]
	_X, _Y = first.(G), last.(G)
	scatter!(_X, _Y, ms=4, color=:white, label=false, marker=:circle)
	if reduced
		xlims!(saved_lims[1]...)
		ylims!(saved_lims[2]...)
	end
	return plot!()
	# return scatter!([p_curr], [v_curr], label=false, color=:white)
end

# ╔═╡ 00874691-d3e0-424e-b0a6-e3457ddcc3cc
md"""
# Utilities
"""

# ╔═╡ b1ea10d1-67be-4c46-a675-519579ab331b
function filltomatch!(A::Vector{Vector{T}}) where T <: Any
	max_length = maximum(map(length, A))	
	for i in 1:length(A)
		len = length(A[i])
		if len < max_length
			# fill last element to `max_length`
			A[i] = vcat(A[i], fill(A[i][end], max_length-len))
		end
	end
	return A
end

# ╔═╡ b1eea5ef-9485-4a55-a213-fb72cf4f2997
function mismatch_mean(A)
	max_length = maximum(map(length, A))
	Z = [map(a->i <= length(a) ? a[i] : nothing, A) for i in 1:max_length]
	return map(mean, map(z->filter(!isnothing, z), Z))
end

# ╔═╡ ca4d62d0-0842-4990-b86b-3e9622d6dee7
function mismatch_std(A)
	max_length = maximum(map(length, A))
	Z = [map(a->i <= length(a) ? a[i] : nothing, A) for i in 1:max_length]
	stds = map(std, map(z->filter(!isnothing, z), Z))
	return map(σ->isnan(σ) ? 0 : σ, stds)
end

# ╔═╡ adb487a5-5ad4-4053-841b-46aeeb2acd5e
function plot_relative_errors(rel_errs_vector; c=:crimson, label=false, hold=false, marker=:circle)
	μ_rel_errs = mismatch_mean(rel_errs_vector)
	σ_rel_errs = mismatch_std(rel_errs_vector)
	ts = [1:length(μ_rel_errs);] .- 1
	plot_fn = hold ? plot! : plot
	plot_fn(ts, μ_rel_errs,
		    title="relative volume error (aggregate: $gds, 10 seeds)",
		    ribbon=σ_rel_errs, fillalpha=0.1,
			marker=(marker, 3, 1.0),
			bg_inside=:white,
		    xlabel="time step", ylabel="relative error", label=label, lw=2, c=c)
	plot!(ts, μ_rel_errs .+ σ_rel_errs, c=c, label=false, alpha=0.33)
	plot!(ts, μ_rel_errs .- σ_rel_errs, c=c, label=false, alpha=0.33)
	annotate!(0.2, 10, text("overestimated", :gray, :left, 10, "Computer Modern"))
	annotate!(0.2, -10, text("underestimated", :gray, :left, 10, "Computer Modern"))
	xlims!(0, Inf)
	# ylims!(minimum(μ_rel_errs), Inf)
	plot!(legend=:bottomright)
end

# ╔═╡ 93545762-a44d-41a8-8e00-e1493d0748e4
begin
	plot_relative_errors(rel_errs_blob;
		label=string("blob: ", accuracy(presults_blob)),
		hold=false, c=colors[1], marker=:star6)
	plot_relative_errors(rel_errs_ellipse;
		label=string("ellipse: ", accuracy(presults_ellipse)),
		hold=true, c=colors[2], marker=:diamond)
	plot_relative_errors(rel_errs_circle;
		label=string("circle: ", accuracy(presults_circle)),
		hold=true, c=colors[3], marker=:circle)
	draw_zero()

	if false
		xlims!(15, 20)
		ylims!(-35, 5)
		title!("")
		# xlabel!("")
		# ylabel!("")
		plot!(legend=false, size=(300,250))
	else
		plot!()
	end
	ylims!(-140, 140)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AverageShiftedHistograms = "77b51b56-6f8f-5c3a-9cb4-d71f9594ea6e"
BSON = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Measures = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
AverageShiftedHistograms = "~0.8.7"
BSON = "~0.3.4"
ColorSchemes = "~3.16.0"
Distributions = "~0.25.46"
Measures = "~0.3.1"
Plots = "~1.25.7"
PlutoUI = "~0.7.33"
StatsBase = "~0.33.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AverageShiftedHistograms]]
deps = ["LinearAlgebra", "RecipesBase", "Statistics", "StatsBase", "UnicodePlots"]
git-tree-sha1 = "8bdad2055f64dd71a25826d752e0222726f25f20"
uuid = "77b51b56-6f8f-5c3a-9cb4-d71f9594ea6e"
version = "0.8.7"

[[deps.BSON]]
git-tree-sha1 = "ebcd6e22d69f21249b7b8668351ebf42d6dc87a1"
uuid = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
version = "0.3.4"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

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

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "6b6f04f93710c71550ec7e16b650c1b9a612d0b6"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.16.0"

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

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "2e97190dfd4382499a4ac349e8d316491c9db341"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.46"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "4a740db447aae0fbeb3ee730de1afbb14ac798a1"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.63.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "aa22e1ee9e722f1da183eb33370df4c1aeb6c2cd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.63.1+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

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

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

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

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

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

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

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

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0b5cfbb704034b5b4c1869e36634438a047df065"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "7e4920a7d4323b8ffc3db184580598450bde8a8e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "da2314d0b0cb518906ea32a497bb4605451811a4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.33"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "37c1631cb3cc36a535105e6d5557864c82cd8c2b"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

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

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e6bf188613555c78062842777b116905a9f9dd49"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.0"

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

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f35e1879a71cca95f4826a14cdbf0b9e253ed918"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.15"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

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

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["Contour", "Crayons", "Dates", "SparseArrays", "StatsBase"]
git-tree-sha1 = "62595983da672758a96f89e07f7fd3735f16c18c"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "2.7.0"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

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

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

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

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

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

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─0239ffc0-3b91-4441-86ae-cf505742d810
# ╠═9bb8e120-8469-11ec-11ba-6d196236bb2f
# ╠═fb741838-9f6d-4920-bc59-c89e01d20d99
# ╠═9af435d0-2301-411f-84fc-cc5b151177e1
# ╠═7ac110de-5ea6-46c9-8a3a-24289b737944
# ╠═1b6df2b6-6536-4bec-aa78-39dfb90c7ad6
# ╟─a67c7a80-8db6-4d6d-b819-1ee46ebc0396
# ╠═5aa567b3-04fa-4aa5-a690-4fb93ccda580
# ╠═095f3331-54e7-472f-aea9-ab32fe45e680
# ╟─602959ed-946d-48b5-b7d8-ef56e644441c
# ╠═f45b80e5-7952-497d-a2af-636ceda95aab
# ╠═df459d63-f23f-48d7-8255-1b3da0910eda
# ╠═70aad8be-1d90-4d01-bfcd-e98206ce89b0
# ╟─2c82b274-f1ba-4368-b990-03407fe58e02
# ╠═3adf7896-4b2b-4421-aee0-0d19eb64f103
# ╠═ea904ef3-03fb-4420-b297-c31372624ed4
# ╠═d3440665-b3e5-453b-9eb7-41a6859ff4b7
# ╠═c7a3e631-16bc-4517-bf82-60e15baaa419
# ╟─d546a273-f8d5-4601-8a1f-4cbbb175ef77
# ╠═a3b8e58a-660a-433a-9dfb-2c75d8de1635
# ╠═3683782b-a884-44de-8986-66c6d7ba9eba
# ╠═182672c0-f623-4a4c-80f4-880afaa69b02
# ╠═4296bde1-e654-481c-b365-61d0ddc6dac7
# ╟─3cbac42c-5c87-41e6-99d5-60eeec2530a5
# ╠═cd5ea37f-5698-47a8-bb79-82f7c65f455d
# ╠═c8fd071a-7f74-45c2-aab0-3fd4c0d65dc8
# ╠═c3352c24-c8fa-44a4-ab22-3e0edee501bf
# ╠═61402b63-a2d7-4d56-8545-aa26ac1a83b3
# ╟─3085e225-59b9-4f0c-9586-5d82f307bf34
# ╠═4ec5cb70-6126-479f-b289-7a79fadcaaf5
# ╠═61b8d873-97f9-40d8-8884-78f0b0f03c7b
# ╠═ecd859be-e8db-4380-a6fd-b85438c63781
# ╠═1192839e-af4e-4ba3-9e3b-a6ebc2fedff6
# ╠═8389a09e-bb8b-4946-854b-0f92b58977a8
# ╠═f0f3efb8-e52d-4552-9675-762c4f3fedc4
# ╟─4e4f6167-875a-41a7-b84a-53d65f5c00a7
# ╠═763aa97a-bbca-4d89-a5f2-843af83baeb8
# ╠═15a60749-135c-4517-af67-3da1cdc62e42
# ╠═6ca61949-4741-48bb-808a-bf41e6eb1c00
# ╠═7f1e37b7-bdbe-41fe-885b-8abcc0e35578
# ╠═c2f7cae6-27ec-4e69-9b3b-410bd63a3443
# ╠═9a9390c0-4946-40a5-8ed0-a85267ba25f9
# ╠═602249bc-4d83-45e3-bd12-cb4097163678
# ╠═9cbb95db-82ec-4053-9624-3ba34d5cd3e1
# ╠═93545762-a44d-41a8-8e00-e1493d0748e4
# ╠═adb487a5-5ad4-4053-841b-46aeeb2acd5e
# ╠═2e4bc9ef-83c7-45ed-a59e-fec98993e53a
# ╠═d972f53b-90fc-4d12-8697-1b18db694a62
# ╟─45ad5118-df2e-44b9-9a24-20067944b06c
# ╠═67db98f5-b5a6-45bd-9676-7b29f8e6f846
# ╠═bcfc1720-5737-4145-aec9-cca110544033
# ╠═dab7bfd4-6e8d-4b88-963b-82647f2cd86e
# ╟─ce02b5d0-4f40-4dcc-897d-8f49082725af
# ╠═499fdc09-fdd6-4ea7-967a-44c4a92a6e4c
# ╠═dfa9ae81-5f0a-41d5-8887-712c11a969a4
# ╠═43b5f33b-120b-4abe-a970-a21fcf697a7f
# ╠═686ebcad-0712-4a2a-ac8f-90389fc5f7b0
# ╠═246dd86d-95bd-42a5-abdb-a7b1280ed42f
# ╠═f5b47f3e-fa1e-4c33-9947-14b8365e293e
# ╠═29de6acb-cb1a-42d7-ae8d-cfeb09be2136
# ╠═f68f7701-746f-4ffa-9a80-6b6d12f85559
# ╠═869bd9fa-3c23-4d0a-b7fa-7e5942b0a556
# ╠═6c5bd780-30f1-404c-ac6c-65a9f775636c
# ╠═f981ef50-efa6-4716-85e8-ef0a90508d79
# ╠═9452720c-95d8-43a5-a399-04589b5d5bde
# ╠═0a458ee4-6781-4f79-94c4-3d5d0d698e9b
# ╠═6c495b97-74b6-4dd5-a66d-2b626899ea2b
# ╠═f0411d8c-a907-4ca7-9d49-a965c884aaff
# ╠═598875f9-e08f-440c-907b-8b840080988c
# ╠═3318952d-9a71-4f7c-bf84-7e3bbf9e1983
# ╠═1b9633f7-f2f6-4de8-9af9-cc325da589fb
# ╟─cf93c16a-5547-4815-9adf-ca0147a445e4
# ╠═0a4cab3f-131f-4741-80b3-5046390581f6
# ╠═7229656e-3525-4be4-a11f-dd277c3b46cd
# ╠═dba66cd7-c91d-4619-8ccc-5f39327ed160
# ╠═02c3f953-69b7-4419-96d2-18341d8f47f8
# ╟─56f9d9a6-82b9-40d4-ba56-fe4768d6b471
# ╠═0d57dadb-0681-4fd9-b834-81bfe0042d47
# ╠═1bf36d36-e6a8-4aca-8d75-0862d4029acc
# ╠═8951fb50-a5b3-4dce-a5d1-2d79f230f28b
# ╠═36476ab4-6c83-49a2-801e-fa442b29f0e1
# ╠═fd340989-63ac-4393-ac7b-da9a0d7de25a
# ╠═c61ac63d-969f-4bcb-a0f3-abbafa764923
# ╠═3f4b47b6-77b2-425b-8ce7-fa38f4bde61f
# ╠═9812856d-eeba-48da-86b4-b8005a577d5e
# ╠═94f09ea2-792a-4936-8f40-292d4fbd376b
# ╟─e9854a55-31ad-4d7f-8331-e5575a6375c6
# ╟─8333262e-8968-4396-8dbf-7a30c5d2f466
# ╠═1b514edc-cde4-4b1f-add8-b21ef545b494
# ╠═5002c00d-c831-4d9a-8c4a-5aaa16ac8819
# ╟─f588a53c-3b72-4f98-89ac-219555d6791c
# ╠═eb14649a-b273-469c-a720-02480641283a
# ╠═7acb3d44-0a71-4a5f-9678-efc716f530d0
# ╟─f85ead75-cc42-4ed5-a410-7c51f6e8a6a5
# ╠═462df3d0-860f-4629-af3e-95648f56c1df
# ╠═b10e3cb2-f8f1-4743-984c-4dc167aa7f59
# ╟─e2052af4-26d8-45ed-aeed-0ab0dd37b6f4
# ╠═d9b4cbe2-884c-407f-843c-f5cbda9c3b14
# ╠═97ff0196-af8e-4060-a590-2784e3a967b1
# ╠═4f04be51-8452-4584-b1c8-804e957de851
# ╠═6570b47d-e84a-4c63-a0d1-a831b647d8d3
# ╟─ce739257-3473-426f-aef4-86a7e016a575
# ╠═cfe230ef-cf77-4a45-9c6b-85781383203d
# ╠═232b0feb-c67d-4ccf-ad2b-bba1802e9344
# ╠═cf562730-8ac6-4b45-a311-a8208c3982fb
# ╠═78d81469-4c41-4be9-afbc-c5ce52e9b2db
# ╠═cf75c166-f776-4b49-b3d2-840aaf8b2231
# ╠═272b238e-c9f6-4856-8d69-1ec8d07c1b0b
# ╠═cc0afc27-926d-4fc7-b802-ea31faaa7041
# ╠═fb8ba111-ce93-4c73-8602-4b6199eba3e0
# ╠═9a02c71a-b8cf-41fc-bca7-490b2c2f8966
# ╠═53e84c72-293a-4818-853a-a9da2e7e580d
# ╠═7a7bc902-5439-4248-9fc3-3baaacab4531
# ╠═59b6f645-616a-4694-9c43-915926d9bf74
# ╠═a5c975e4-8118-421c-a55e-a4d438337c71
# ╟─83a79501-d3ba-4adf-b183-4d1b28b815bd
# ╠═a3ba5f18-7d8e-4a57-85e1-56933771bd77
# ╠═fcf97cb8-b9c4-40f6-a694-af848046d4bd
# ╠═94364d2e-7936-4919-9641-5d1ec1762f36
# ╠═dc4e45cc-eb23-4dc5-83dc-a2ab136db691
# ╠═e7968ebc-8a30-4a9c-b260-fc86883de5a9
# ╠═989dea24-8610-4275-b4fb-a601a4bc870d
# ╠═1e3318e5-3181-4d8a-8602-22291b3750ef
# ╟─53d2a31e-0fd5-49a9-ae81-5a55b22e328d
# ╠═4d262ca1-e5e9-4009-a445-59395dd25503
# ╠═ac753e4d-e7f6-4eb6-af22-d06f58473fb7
# ╠═d1188bf0-92c1-4790-a6a6-953a18f46354
# ╠═abd7e135-f4b5-4cbe-86da-19e94039a7e7
# ╠═c623ed18-c95c-4ef3-88c8-57a1213660bf
# ╠═d76d272b-d810-4726-86ab-fa89f81bcb12
# ╠═f57dde29-d588-4fa5-bd41-f9b13b9d877e
# ╠═4a552377-db94-47df-a438-2309bcb7b814
# ╠═74fdc433-91ff-4c08-9c46-2563be98f317
# ╠═eaebf4df-c915-4136-aa99-7101563d4aca
# ╠═b3296ce7-f606-4c43-b10a-3603c887be4e
# ╠═94b106f6-409d-457f-a4a7-4ad68a9fd0ef
# ╠═3bad29de-edc2-4248-b5fb-5c1aa360c457
# ╠═c1bccfa7-9f45-4ed2-8a4b-5a269329ed0a
# ╠═ea29156a-4e96-4f7c-a9ca-af235259d14b
# ╠═04144789-7339-4f3e-943b-957d87e6dc33
# ╠═7017212c-71da-4fe2-90c6-fe95534446da
# ╟─c21fe459-c5c6-4324-ac17-005dddc94bb7
# ╠═d06f41df-45ee-4ed0-adb1-ddd5666bd243
# ╠═d79d7b50-6714-4175-8eba-ef3ad7451de9
# ╠═c69f4e20-0489-4d93-98e0-e47a78a9fb55
# ╠═2563b2e6-d2c6-4386-880d-48e8ba31cd44
# ╠═e3d813b2-1c52-4d05-8953-50ed8412cc6e
# ╠═bf75df30-32cf-49c0-8048-b891f5f03e13
# ╠═94a167e5-7a24-4acc-8fdf-3fc5ad3f4af7
# ╟─b7d46821-abb0-42ba-88a1-f9e24371a455
# ╠═207d023f-bd15-4e73-9de8-6291dea3ec39
# ╠═600eae30-53bf-46b1-863e-ca28b64d1834
# ╠═4f14d88e-7ef1-4ab2-992b-5afbaddb3e37
# ╠═367369f3-c2cb-4c56-b4fe-e29c9536289f
# ╠═39803ea6-ecd5-48c4-a4b7-729ab2743263
# ╠═4d1911a1-2c79-4524-a630-5eaf9e462ffe
# ╠═c57fddaa-5e69-4ea0-8469-1784ab1019bb
# ╠═0b91ab7d-f5bf-4889-a7ee-d1b5783e72e2
# ╠═914a81ee-e8c5-4173-bceb-710dbfedc4a3
# ╠═08102844-7ba0-4ca2-97fc-7b66410859f4
# ╠═d8723841-f4b2-4c09-bac3-3295b8e5176d
# ╠═166ad666-5a03-410c-9448-30c518051d3b
# ╟─7b56bf68-6d3c-4b6e-8826-082215409208
# ╠═0e39f099-62d3-49ba-814d-9df55be91f81
# ╠═ba97c361-013a-41f3-a54c-02460ef63270
# ╠═10e845f2-2551-4453-a7c2-180118503283
# ╠═172f1589-5a7c-4dd9-aac4-3be044dd32df
# ╠═73abd22c-7438-4398-bc43-6441ca83939f
# ╠═45da8856-59b3-4dab-9d9f-edae3269baad
# ╠═70b883ec-15b7-432c-b375-d5519ea32028
# ╠═5c84a2ac-b719-47fd-b406-81a2255aa9c2
# ╠═0ad9c8b0-918f-4e4c-b9a9-b41e296d7dd4
# ╠═9d7caa02-54a1-40f0-816a-85faff5992d6
# ╠═532beb9a-06c0-4ef7-b905-813aaec5a55d
# ╠═ef812da1-135c-4664-b79b-91430f46d3dc
# ╠═9e169ad0-b323-4cd8-a60f-fbbd35bca36a
# ╠═3afb0730-51cb-4f7f-a26b-8618df4126d8
# ╠═19236596-7225-4767-8348-4be4caa684c2
# ╠═5ad33171-5466-4165-9615-6b618eb4f5f8
# ╠═96182d6b-4821-4abc-8a41-3be25b9f6e0b
# ╟─00874691-d3e0-424e-b0a6-e3457ddcc3cc
# ╠═b1ea10d1-67be-4c46-a675-519579ab331b
# ╠═b1eea5ef-9485-4a55-a213-fb72cf4f2997
# ╠═ca4d62d0-0842-4990-b86b-3e9622d6dee7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
