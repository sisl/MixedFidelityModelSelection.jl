### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ e4de635f-3564-4044-a62d-d7ee9a6ba95c
begin
	using Revise
	using Pkg
	Pkg.develop(path="..//")
	using MixedFidelityModelSelection
	using BSON
	using Statistics
	using PlutoUI
	using Plots; default(fontfamily="Computer Modern", framestyle=:box)
end

# ╔═╡ 763aa97a-bbca-4d89-a5f2-843af83baeb8
using ColorSchemes

# ╔═╡ e7968ebc-8a30-4a9c-b260-fc86883de5a9
using Measures

# ╔═╡ 4d262ca1-e5e9-4009-a445-59395dd25503
using StatsBase

# ╔═╡ 2563b2e6-d2c6-4386-880d-48e8ba31cd44
using LinearAlgebra

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

# ╔═╡ 2f0b58fb-5857-4dfe-b1e7-bf4a042e78cb
md"""
# New data format
"""

# ╔═╡ ec1e38f4-6e90-4041-8a96-b1d829de193c
results = BSON.load("..\\scripts\\results\\results_fixed_bores.bson")[:results]

# ╔═╡ 9085cf37-0390-482f-94b4-40e46ce3d51e
begin
	results50x50blob = results[(:blob, (50,50,1), 100)]
	results50x50ellipse = results[(:ellipse, (50,50,1), 100)]
	results50x50circle = results[(:circle, (50,50,1), 100)]
    results30x30blob = results[(:blob, (30,30,1), 100)]
    results30x30ellipse = results[(:ellipse, (30,30,1), 100)]
    results30x30circle = results[(:circle, (30,30,1), 100)]
end;

# ╔═╡ 13fc2600-3789-4fc7-a31e-357f35db4b37
begin
	ptiming_blob = map(t->t.time, results50x50blob[:timing])
	ptiming_ellipse = map(t->t.time, results50x50ellipse[:timing])
	ptiming_circle = map(t->t.time, results50x50circle[:timing])
end;

# ╔═╡ 2c82b274-f1ba-4368-b990-03407fe58e02
md"""
## Loading data (blob)
"""

# ╔═╡ 3adf7896-4b2b-4421-aee0-0d19eb64f103
analysis_blob = BSON.load("..\\results\\$subdir\\analysis_blob_$(gds)_pomcp.100.bson")[:analysis];

# ╔═╡ ea904ef3-03fb-4420-b297-c31372624ed4
presults_blob = analysis_blob["results"]

# ╔═╡ d546a273-f8d5-4601-8a1f-4cbbb175ef77
md"""
## Loading data (ellipse)
"""

# ╔═╡ a3b8e58a-660a-433a-9dfb-2c75d8de1635
analysis_ellipse = BSON.load("..\\results\\$subdir\\analysis_ellipse_$(gds)_pomcp.100.bson")[:analysis];

# ╔═╡ 3683782b-a884-44de-8986-66c6d7ba9eba
presults_ellipse = analysis_ellipse["results"]

# ╔═╡ 3cbac42c-5c87-41e6-99d5-60eeec2530a5
md"""
## Loading data (circle)
"""

# ╔═╡ cd5ea37f-5698-47a8-bb79-82f7c65f455d
analysis_circle = BSON.load("..\\results\\$subdir\\analysis_circle_$(gds)_pomcp.100.bson")[:analysis];

# ╔═╡ c8fd071a-7f74-45c2-aab0-3fd4c0d65dc8
presults_circle = analysis_circle["results"]

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

# ╔═╡ b8a2f550-43f8-4379-8cec-a0c4522c2ab9
md"""
## Colors
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

# ╔═╡ 4e4f6167-875a-41a7-b84a-53d65f5c00a7
md"""
## Relative volume errors
"""

# ╔═╡ 6f2c7efc-81f4-4d79-a58b-5761190b3e7b
MixedFidelityModelSelection.annotate!

# ╔═╡ 93545762-a44d-41a8-8e00-e1493d0748e4
plot_rel_error_aggregate(results50x50blob, results50x50ellipse, results50x50circle)

# ╔═╡ 378681c1-8754-4381-bf94-fa956dc58e70
plot_rel_error_aggregate(results30x30blob, results30x30ellipse, results30x30circle)

# ╔═╡ 45ad5118-df2e-44b9-9a24-20067944b06c
md"""
## Polar comparison plots
- runtime
- regret
- true positive
- false positive
"""

# ╔═╡ ce02b5d0-4f40-4dcc-897d-8f49082725af
md"""
## Confusion statistics
"""

# ╔═╡ 7229656e-3525-4be4-a11f-dd277c3b46cd
cm = confusion_matrix(decisions, true_decisions)

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
regrets = regret(data5[1], 10, 10)

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
	plot(x->ecdf(regrets)(x), label=false, fill=true, 0, 150, size=(300,200), fillcolor=:gray, c=:black, title="cumulative regret (10×10, 10 iters.)", titlefontsize =10)
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

# ╔═╡ Cell order:
# ╟─0239ffc0-3b91-4441-86ae-cf505742d810
# ╠═e4de635f-3564-4044-a62d-d7ee9a6ba95c
# ╠═1b6df2b6-6536-4bec-aa78-39dfb90c7ad6
# ╟─a67c7a80-8db6-4d6d-b819-1ee46ebc0396
# ╠═5aa567b3-04fa-4aa5-a690-4fb93ccda580
# ╠═095f3331-54e7-472f-aea9-ab32fe45e680
# ╟─602959ed-946d-48b5-b7d8-ef56e644441c
# ╠═f45b80e5-7952-497d-a2af-636ceda95aab
# ╠═df459d63-f23f-48d7-8255-1b3da0910eda
# ╠═70aad8be-1d90-4d01-bfcd-e98206ce89b0
# ╟─2f0b58fb-5857-4dfe-b1e7-bf4a042e78cb
# ╠═ec1e38f4-6e90-4041-8a96-b1d829de193c
# ╠═9085cf37-0390-482f-94b4-40e46ce3d51e
# ╠═13fc2600-3789-4fc7-a31e-357f35db4b37
# ╟─2c82b274-f1ba-4368-b990-03407fe58e02
# ╠═3adf7896-4b2b-4421-aee0-0d19eb64f103
# ╠═ea904ef3-03fb-4420-b297-c31372624ed4
# ╟─d546a273-f8d5-4601-8a1f-4cbbb175ef77
# ╠═a3b8e58a-660a-433a-9dfb-2c75d8de1635
# ╠═3683782b-a884-44de-8986-66c6d7ba9eba
# ╟─3cbac42c-5c87-41e6-99d5-60eeec2530a5
# ╠═cd5ea37f-5698-47a8-bb79-82f7c65f455d
# ╠═c8fd071a-7f74-45c2-aab0-3fd4c0d65dc8
# ╟─3085e225-59b9-4f0c-9586-5d82f307bf34
# ╠═4ec5cb70-6126-479f-b289-7a79fadcaaf5
# ╠═61b8d873-97f9-40d8-8884-78f0b0f03c7b
# ╠═ecd859be-e8db-4380-a6fd-b85438c63781
# ╠═1192839e-af4e-4ba3-9e3b-a6ebc2fedff6
# ╠═8389a09e-bb8b-4946-854b-0f92b58977a8
# ╠═f0f3efb8-e52d-4552-9675-762c4f3fedc4
# ╟─b8a2f550-43f8-4379-8cec-a0c4522c2ab9
# ╠═763aa97a-bbca-4d89-a5f2-843af83baeb8
# ╠═15a60749-135c-4517-af67-3da1cdc62e42
# ╠═6ca61949-4741-48bb-808a-bf41e6eb1c00
# ╠═7f1e37b7-bdbe-41fe-885b-8abcc0e35578
# ╟─4e4f6167-875a-41a7-b84a-53d65f5c00a7
# ╠═6f2c7efc-81f4-4d79-a58b-5761190b3e7b
# ╠═93545762-a44d-41a8-8e00-e1493d0748e4
# ╠═378681c1-8754-4381-bf94-fa956dc58e70
# ╟─45ad5118-df2e-44b9-9a24-20067944b06c
# ╟─ce02b5d0-4f40-4dcc-897d-8f49082725af
# ╠═7229656e-3525-4be4-a11f-dd277c3b46cd
# ╠═02c3f953-69b7-4419-96d2-18341d8f47f8
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
