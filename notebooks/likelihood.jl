### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 41dc8c80-9420-11ec-1606-a52cb5907d71
begin
	using Revise
	using Pkg

	Pkg.develop(path="..//")
	using POMDPModelFidelityFramework
	Pkg.develop(path="..//..//MineralExploration//")
	using MineralExploration
	Pkg.develop(path="..//scripts//MEParallel.jl//")
	using MEParallel

	using BSON
	using Statistics
	using PlutoUI
	using Plots; default(fontfamily="Computer Modern", framestyle=:box)
end

# ╔═╡ fb539315-3f2b-40f1-a75f-d470dba2bb93
using Distributions

# ╔═╡ 34ce6bb9-db0a-4850-9bdb-55beb88bf2cf
using LinearAlgebra

# ╔═╡ 6fc1ea28-b663-4ace-8b1f-1d2c574aef2d
using POMDPSimulators

# ╔═╡ 620b9753-e132-49fa-8a82-18468d9fa6df
using LaTeXStrings

# ╔═╡ 16c5b0f2-a15d-4829-b924-f14802b726fd
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 9663e639-14ee-49c8-87c3-4619d2b03db3
config = MEConfiguration(1, (30,30,1), 10, CircleNode, MEJobParameters(name="test", min_bores=15, max_bores=15))

# ╔═╡ 098fb798-05e4-4c3c-bf79-a1822f4c556d
trial = MEParallel.initialize(config)

# ╔═╡ b3aea315-3af2-41de-b767-03f85fbf1e81
b0 = trial.b0;

# ╔═╡ adac57f4-dd06-404c-8e56-91bf484bd9bb
MineralExploration.plot(b0)

# ╔═╡ d90540e9-a8de-4002-8d56-ed136d81466b
a = MineralExploration.MEAction(:drill, CartesianIndex(10,10))

# ╔═╡ b39af824-e482-4cd1-a52f-88defd9efec0
(sp, o, r) = MineralExploration.@gen(:sp,:o,:r)(trial.problem, trial.s0, a, MineralExploration.Random.GLOBAL_RNG)

# ╔═╡ f7eb845b-c8ee-4fb0-98de-3d9f9aef6f6a
bp = MineralExploration.update(trial.updater, b0, a, o)

# ╔═╡ 8858f516-c408-4ed9-974e-83eba4ee5d4b
MineralExploration.plot(bp)

# ╔═╡ aecdb341-635b-4d85-8158-0b5f2ac841d5
bp.particles[1].ore_map |> MineralExploration.plot_ore_map

# ╔═╡ 144e5693-9624-413d-b23e-9736b7312a28
bp.particles[1].ore_map

# ╔═╡ f70064a6-8c71-4354-8203-64773582294e
bp

# ╔═╡ 0dbb781c-5c28-4beb-a795-ad89105db513
MineralExploration.summarize(bp)

# ╔═╡ fad58115-e5a7-480a-b4a5-67d398597475
b0.rock_obs

# ╔═╡ dd9fbe71-f1b7-4dd6-b2a2-1259870fc6d8
md"""
$p(s \mid b_0)$
"""

# ╔═╡ 6ca24aad-28d8-41a8-9b02-7d33810f2f5e
MineralExploration.summarize(b0)[1] |> MineralExploration.plot_ore_map

# ╔═╡ 187fccf2-57f4-41c6-81c6-09506dce5612
function summarize_mainbody(b::MEBelief)
    (x, y, z) = size(b.particles[1].mainbody_map)
    μ = zeros(Float64, x, y, z)
    w = 1.0/length(b.particles)
    for p in b.particles
        mainbody_map = p.mainbody_map
        μ .+= mainbody_map .* w
    end
    σ² = zeros(Float64, x, y, z)
    for p in b.particles
        mainbody_map = p.mainbody_map
        σ² .+= w*(mainbody_map - μ).^2
    end
    return (μ, σ²)
end

# ╔═╡ 869590a3-f217-4bfb-839a-3b13a2131276
summarize_mainbody(b0)[1] |> MineralExploration.plot_ore_map

# ╔═╡ 40491708-07ac-4115-8cfe-ca0e21bf22b9
b0.particles[1].mainbody_map |> MineralExploration.plot_ore_map

# ╔═╡ 09d331e6-3af3-41b0-91f8-ed3ebfcc6c93
function likelihood(b::MEBelief, prior::MEBelief)
	ws = Float64[]
	rock_obs = b.rock_obs
    bore_coords = rock_obs.coordinates
    n = size(bore_coords)[2]
    ore_obs = [o for o in rock_obs.ore_quals]
    K = MineralExploration.calc_K(b.geostats, rock_obs)
    mu = zeros(Float64, n) .+ trial.updater.m.gp_mean
    gp_dist = MvNormal(mu, K)
    # for s in b.particles
        for sᵦ in prior.particles
	        o_n = zeros(Float64, n)
			for i = 1:n
				o_mainbody = sᵦ.mainbody_map[bore_coords[1, i], bore_coords[2, i]]
				o_n[i] = ore_obs[i] - o_mainbody
			end
			w = pdf(gp_dist, o_n)
			# TODO: likelihood of mainbody params (change how belief is sampled)
			push!(ws, w)
        end
    # end
    return ws
end

# ╔═╡ cdc1dd82-2106-4d7c-9125-a62c7f031b9a
mainbody_gen = CircleNode(grid_dims=(30,30,1))

# ╔═╡ 5d5fa2ea-33a6-4241-9e54-b104b027da55
mainbody_map, mainbody_params = rand(mainbody_gen)

# ╔═╡ d1c1dd09-cd85-41ce-b7e0-52300191214c
MineralExploration.plot_ore_map(mainbody_map)

# ╔═╡ e3368e6c-65d3-465b-8ad0-101d075e52e1
mainbody_params

# ╔═╡ e59bc42c-e65f-4e60-ba75-dd6c7d2cf9a4
function llhood(mainbody_gen::CircleNode, params) # TODO: likelihood
	center, radius, _ = params
	pdf_center = pdf(mainbody_gen.center, center)
	pdf_radius = pdf(mainbody_gen.radius, radius)
	return pdf_center*pdf_radius
end

# ╔═╡ e3ee4a87-2219-48de-8e74-5328fdb08557
llhood(mainbody_gen, mainbody_params)

# ╔═╡ 692dd95f-18f9-4c95-b10b-d3df4e9ff6c8
function likelihood2(b::MEBelief, prior::MEBelief)
	ws = Float64[]
	rock_obs = b.rock_obs
    bore_coords = rock_obs.coordinates
    n = size(bore_coords, 2)
    ore_obs = [o for o in rock_obs.ore_quals]
    K = MineralExploration.calc_K(b.geostats, rock_obs)
    mu = zeros(Float64, n) .+ trial.updater.m.gp_mean
    gp_dist = MvNormal(mu, K)
	for sᵦ in prior.particles
		o_n = zeros(Float64, n)
		for i = 1:n
			o_mainbody = sᵦ.mainbody_map[bore_coords[1, i], bore_coords[2, i]]
			o_n[i] = ore_obs[i] - o_mainbody
		end
		w = pdf(gp_dist, o_n) # * llhood(mainbody_gen, mainbody_params)
		# TODO: likelihood of mainbody params (change how belief is sampled)
		push!(ws, w)
	end
    return ws
end

# ╔═╡ 90912736-05d2-402e-ba41-486e658128e0
ℓ = likelihood2(bp, b0)

# ╔═╡ 6efa83e1-b49a-4e60-93f7-a8e21b746a71
length(ℓ)

# ╔═╡ f4d561df-f294-4a02-8de9-f80bcf39dbe5
mean(ℓ)

# ╔═╡ c8e345eb-1afa-4010-8155-aeba477827f8
histogram(ℓ)

# ╔═╡ f82d7626-2495-4779-ad26-57291e45e5dd
arg = argmax(likelihood(bp, b0))

# ╔═╡ e7b3e4a4-56c8-42d3-9968-54ddf45cddbc
md"""
# Truth vs. belief
"""

# ╔═╡ 33ff7adc-f6d3-4839-92ed-520229366a7a
if true
	L = Float64[]
	for (sp,a,o,bp,t) in stepthrough(trial.problem, trial.planner, trial.updater, b0, trial.s0, "sp,a,o,bp,t", max_steps=50)
		@info t
		ℓ = mean(likelihood2(bp, b0))
		push!(L, ℓ)
	end
	plot(log.(L), legend=:topleft)
end

# ╔═╡ 52fa0795-34a5-4fb2-9f00-a9ed8e8203af
plot(L)

# ╔═╡ 988e0df8-cba9-439f-8188-9e42752ad277
trial.s0.ore_map |> MineralExploration.plot_ore_map

# ╔═╡ 17b76c3c-5242-412a-83d0-128cea37e39b
trial.s0.ore_map[a.coords[1], a.coords[2], 1]

# ╔═╡ baabd323-4340-4f79-9d5b-26701216d9a8
begin
	bp.particles[arg].ore_map |> MineralExploration.plot_ore_map
	title!("belief particle")
end

# ╔═╡ 8dde997d-85bb-4d42-8c62-308a6737c40e
bp.particles[arg].ore_map[a.coords[1], a.coords[2], 1]

# ╔═╡ 1594de68-316e-469f-b5d3-6154f14b421a
MEParallel.initialize(MEConfiguration(1, (30,30,1), 10, BlobNode,
			MEJobParameters(name="test", min_bores=10, max_bores=10))).s0.ore_map |> MineralExploration.plot_ore_map

# ╔═╡ 2c70f2e0-7201-401a-9731-ac7f509d9caf
if true
	μ_L = Dict()
	σ_L = Dict()
	b_final = Dict()
	for shape in [BlobNode, EllipseNode, CircleNode]
		config = MEConfiguration(1, (10,10,1), 10, shape,
			MEJobParameters(name="test", min_bores=10, max_bores=10))
		trial = MEParallel.initialize(config)
		μ_L[shape] = Float64[]
		σ_L[shape] = Float64[]
		
		for (sp,a,o,bp,t) in stepthrough(trial.problem, trial.planner, trial.updater,
			                             trial.b0, trial.s0, "sp,a,o,bp,t",
			                             max_steps=50)
			@info t
			ℓ = likelihood(bp, trial.b0)
			push!(μ_L[shape], mean(ℓ))
			push!(σ_L[shape], std(ℓ))
			b_final[shape] = bp
		end
	end
	μ_L, σ_L
end

# ╔═╡ bf220a69-2e3a-4c9d-bb8f-4c097d0d8c99
plot(b_final[BlobNode])

# ╔═╡ 917be533-532e-45ce-a4ac-affe13cc5240
plot(b_final[EllipseNode])

# ╔═╡ 4b571bf3-fe5c-42c6-abdd-716bd8435d77
plot(b_final[CircleNode])

# ╔═╡ 1a6acfb4-f108-4a05-a205-af03b3d58a74
begin
	pgfplotsx()
	plot(log.(μ_L[BlobNode]), ribbon=log.(σ_L[BlobNode])/2, fillalpha=0.15,
		label="blob", lw=2, ls=:dash, marker=:x, ms=3) # :star6
	plot!(log.(μ_L[EllipseNode]), ribbon=log.(σ_L[EllipseNode])/2, fillalpha=0.15,
		label="ellipse", lw=2, ls=:dash, marker=:diamond, ms=3)
	plot!(log.(μ_L[CircleNode]), ribbon=log.(σ_L[CircleNode])/2, fillalpha=0.15,
		label="circle", lw=2, ls=:dash, marker=:circle, ms=3)
	ylabel!("log-likelihood")
	xlabel!("time step")
	plot!(legend=:bottomleft, size=(400,200))
	plt = title!(L"\log\; p(b_t \mid b_0)")
	gr()
	plt
end

# ╔═╡ 3e7706fc-a91a-4731-be48-c8cb6e9aa7c6
savefig(plt, "likelihood.tex")

# ╔═╡ edc082cd-c2c4-446d-a265-126ac72b8c25
Plots.supported_markers()

# ╔═╡ 9ae3326c-e242-42b1-b2f1-ec28e3efe881
function bounds(X)
	[minimum(X), maximum(X)]
end

# ╔═╡ 03b3b5b9-1f64-4a77-afaf-c96d67b7ea11
bounds(log.(μ_L[BlobNode])) |> diff

# ╔═╡ 1ebcaf61-f038-4209-aa20-ca9e65cd197c
bounds(log.(μ_L[EllipseNode])) |> diff

# ╔═╡ d1893c10-e5ef-4d10-a8b3-9b1f0501f808
bounds(log.(μ_L[CircleNode])) |> diff

# ╔═╡ 63c886d5-20f2-420e-ad48-77558de04795
begin
	plot(diff(log.(μ_L[BlobNode])), ribbon=diff(log.(σ_L[BlobNode])/2), fillalpha=0.15,	label="blob", lw=2, ls=:dash, marker=:star6, ms=3)
	plot!(diff(log.(μ_L[EllipseNode])), ribbon=diff(log.(σ_L[EllipseNode])/2), fillalpha=0.15,	label="ellipse", lw=2, ls=:dash, marker=:diamond, ms=3)
	plot!(diff(log.(μ_L[CircleNode])), ribbon=diff(log.(σ_L[CircleNode])/2), fillalpha=0.15,
		label="circle", lw=2, ls=:dash, marker=:circle, ms=3)
	ylabel!("log-likelihood gradient")
	xlabel!("time step")
	plot!(legend=:bottomleft)
	title!(L"\nabla\log\; p(b' \mid b_0)")
end

# ╔═╡ 34fdd973-a081-4521-8172-664f1b9641df
diff(log.(μ_L[CircleNode]))

# ╔═╡ 5c810a1a-2f0b-4cc5-9a05-2c5975a74cd7
md"""
---
"""

# ╔═╡ ad371651-c9ff-4ac7-bb4c-a977d8ed2327
function linear_regression(X, y)
    𝐗 = [ones(size(y)) X]
    θ = 𝐗\y
    return x -> [1;x]'θ
end

# ╔═╡ 2a2db048-1581-46ab-890e-af9c4a545d12
X = 1:10

# ╔═╡ 6de07d23-bcbd-4ac8-84b8-dbc818387ffc
y = [X;] + randn(size(X,1));

# ╔═╡ f3679948-0c5e-4f0a-888b-db6aa614f298
f = linear_regression(X,y)

# ╔═╡ 1b78fd90-edd0-4225-b6fe-406ba7a4d030
begin
	plot(X, f, label="fit", color=:blue, legend=:topleft)
	scatter!(X, y, label="data", color=:red, title="linear regression fit")
end

# ╔═╡ 859cd515-0b34-460a-bf08-33e43f743695
θ = [4.67, 5.44]

# ╔═╡ c3650a56-7479-43c0-be7d-4d31d6ce1f06
function get_parameters(f)
	intercept = f(0)
	slope = f(1) - intercept
	return (intercept, slope)
end

# ╔═╡ 95d7c96c-9a8b-40d1-871a-e3468ad113ce
get_parameters(x->[1;x]'θ)

# ╔═╡ cceffc6a-458f-47fa-832a-6f5711a706b1
md"""
# Step through (quickmath)
"""

# ╔═╡ 31091273-f1d4-4a67-935f-2526e0841d98
b0

# ╔═╡ 2e379c61-b85f-4ace-bc1a-6b04c6452b56
if true
	for (sp,a,o,bp,t) in stepthrough(trial.problem, trial.planner, trial.updater, b0, trial.s0, "sp,a,o,bp,t", max_steps=5)
		@info "" plot(bp)
	end
end

# ╔═╡ bec867c7-2227-4a5f-aff1-b0a1bc65b398
md"""
# Belief update
"""

# ╔═╡ 4fdd5a9f-cc78-40a1-bd4b-343e964c4a1f
b′ = MineralExploration.update(trial.updater, trial.b0, a, o);

# ╔═╡ eb42e483-5c95-42e9-9211-a47fedfdadd0
plot(b′)

# ╔═╡ 0142f9a5-a28f-4f70-8207-f1300ba86c5d
a2 = MineralExploration.MEAction(:drill, CartesianIndex(18,18))

# ╔═╡ 6fac34c4-2711-4105-b3df-bf063a166b3b
(sp2, o2, _) = MineralExploration.@gen(:sp,:o,:r)(trial.problem, sp, a2, MineralExploration.Random.GLOBAL_RNG);

# ╔═╡ 3ed4009b-0f82-470c-b8ec-df3183f7fd66
o2.ore_quality

# ╔═╡ d3102841-af94-4f70-a90d-0f21d80d193e
b2′ = MineralExploration.update(trial.updater, b′, a2, o2);

# ╔═╡ c38b74ce-f51a-4508-9f71-15a11135aefc
plot(b2′)

# ╔═╡ 56438479-98e8-4793-a2f9-6bbc555b7993
MineralExploration.GeoStats.sill(b0.geostats.variogram)

# ╔═╡ 8c168dc9-6a16-4baf-a2c7-20d7bdc59772
domain = MineralExploration.PointSet(sp2.rock_obs.coordinates)

# ╔═╡ 15e5185d-a1a7-46e7-a9da-1650e72b3ad0
Dd = [MineralExploration.GeoStats.Point(p.coords[1] + 0.5, p.coords[2] + 0.5) for p in domain]

# ╔═╡ b0bc94b4-f65c-41fb-a72a-bf67cf1c6076
MineralExploration.GeoStats.pairwise(b0.geostats.variogram, Dd)

# ╔═╡ 016eced6-397f-4c68-a0d8-8de03298dffb
K = MineralExploration.GeoStats.sill(b0.geostats.variogram) .-
	MineralExploration.GeoStats.pairwise(b0.geostats.variogram, Dd)

# ╔═╡ d5454c4e-748c-4888-9102-eb6cefbc8409
b0.geostats.variogram

# ╔═╡ 4061e3bf-0b1c-48b5-abf2-76dd3b371f52
begin
	plot(b0.geostats.variogram, lw=2)
	hline!([b0.geostats.variogram.sill], label="sill", color=:red, ls=:dash)
	vline!([b0.geostats.variogram.range], label="range", color=:green, ls=:dash)
	scatter!([0], [b0.geostats.variogram.nugget], label="nugget", color=:gold)
	plot!(legend=:bottomright,
		  ylims=(-0.0002, 0.0055), xlims=(-2, b0.geostats.variogram.range*3))
	hline!([0], color=:gray, label=false)
	vline!([0], color=:gray, label=false)
end

# ╔═╡ 6980e00b-8bb3-4486-94b7-43095a91c4eb
begin
	n = size(sp2.rock_obs.coordinates,2)
	μ = zeros(n) .+ trial.problem.gp_mean
end

# ╔═╡ a490b4b1-c21b-41bb-8915-3b2e215df15c
gp_dist = MvNormal(μ, K)

# ╔═╡ 81c75b21-6910-47fa-988d-8f240d9e1473
begin
	ranges = 0:0.01:1
	heatmap(ranges, ranges, (x,y)->pdf(gp_dist, [x,y]), ratio=1, c=:viridis)
	xlims!(first(ranges), last(ranges))
end

# ╔═╡ Cell order:
# ╠═41dc8c80-9420-11ec-1606-a52cb5907d71
# ╟─16c5b0f2-a15d-4829-b924-f14802b726fd
# ╠═9663e639-14ee-49c8-87c3-4619d2b03db3
# ╠═098fb798-05e4-4c3c-bf79-a1822f4c556d
# ╠═b3aea315-3af2-41de-b767-03f85fbf1e81
# ╠═adac57f4-dd06-404c-8e56-91bf484bd9bb
# ╠═d90540e9-a8de-4002-8d56-ed136d81466b
# ╠═b39af824-e482-4cd1-a52f-88defd9efec0
# ╠═f7eb845b-c8ee-4fb0-98de-3d9f9aef6f6a
# ╠═8858f516-c408-4ed9-974e-83eba4ee5d4b
# ╠═aecdb341-635b-4d85-8158-0b5f2ac841d5
# ╠═144e5693-9624-413d-b23e-9736b7312a28
# ╠═f70064a6-8c71-4354-8203-64773582294e
# ╠═0dbb781c-5c28-4beb-a795-ad89105db513
# ╠═fad58115-e5a7-480a-b4a5-67d398597475
# ╠═fb539315-3f2b-40f1-a75f-d470dba2bb93
# ╠═dd9fbe71-f1b7-4dd6-b2a2-1259870fc6d8
# ╠═6ca24aad-28d8-41a8-9b02-7d33810f2f5e
# ╠═187fccf2-57f4-41c6-81c6-09506dce5612
# ╠═869590a3-f217-4bfb-839a-3b13a2131276
# ╠═40491708-07ac-4115-8cfe-ca0e21bf22b9
# ╠═09d331e6-3af3-41b0-91f8-ed3ebfcc6c93
# ╠═cdc1dd82-2106-4d7c-9125-a62c7f031b9a
# ╠═5d5fa2ea-33a6-4241-9e54-b104b027da55
# ╠═d1c1dd09-cd85-41ce-b7e0-52300191214c
# ╠═e3368e6c-65d3-465b-8ad0-101d075e52e1
# ╠═e59bc42c-e65f-4e60-ba75-dd6c7d2cf9a4
# ╠═e3ee4a87-2219-48de-8e74-5328fdb08557
# ╠═34ce6bb9-db0a-4850-9bdb-55beb88bf2cf
# ╠═692dd95f-18f9-4c95-b10b-d3df4e9ff6c8
# ╠═90912736-05d2-402e-ba41-486e658128e0
# ╠═6efa83e1-b49a-4e60-93f7-a8e21b746a71
# ╠═f4d561df-f294-4a02-8de9-f80bcf39dbe5
# ╠═c8e345eb-1afa-4010-8155-aeba477827f8
# ╠═f82d7626-2495-4779-ad26-57291e45e5dd
# ╟─e7b3e4a4-56c8-42d3-9968-54ddf45cddbc
# ╠═33ff7adc-f6d3-4839-92ed-520229366a7a
# ╠═52fa0795-34a5-4fb2-9f00-a9ed8e8203af
# ╠═988e0df8-cba9-439f-8188-9e42752ad277
# ╠═17b76c3c-5242-412a-83d0-128cea37e39b
# ╠═baabd323-4340-4f79-9d5b-26701216d9a8
# ╠═8dde997d-85bb-4d42-8c62-308a6737c40e
# ╠═6fc1ea28-b663-4ace-8b1f-1d2c574aef2d
# ╠═1594de68-316e-469f-b5d3-6154f14b421a
# ╠═2c70f2e0-7201-401a-9731-ac7f509d9caf
# ╠═bf220a69-2e3a-4c9d-bb8f-4c097d0d8c99
# ╠═917be533-532e-45ce-a4ac-affe13cc5240
# ╠═4b571bf3-fe5c-42c6-abdd-716bd8435d77
# ╠═620b9753-e132-49fa-8a82-18468d9fa6df
# ╠═1a6acfb4-f108-4a05-a205-af03b3d58a74
# ╠═3e7706fc-a91a-4731-be48-c8cb6e9aa7c6
# ╠═edc082cd-c2c4-446d-a265-126ac72b8c25
# ╠═9ae3326c-e242-42b1-b2f1-ec28e3efe881
# ╠═03b3b5b9-1f64-4a77-afaf-c96d67b7ea11
# ╠═1ebcaf61-f038-4209-aa20-ca9e65cd197c
# ╠═d1893c10-e5ef-4d10-a8b3-9b1f0501f808
# ╠═63c886d5-20f2-420e-ad48-77558de04795
# ╠═34fdd973-a081-4521-8172-664f1b9641df
# ╟─5c810a1a-2f0b-4cc5-9a05-2c5975a74cd7
# ╠═ad371651-c9ff-4ac7-bb4c-a977d8ed2327
# ╠═2a2db048-1581-46ab-890e-af9c4a545d12
# ╠═6de07d23-bcbd-4ac8-84b8-dbc818387ffc
# ╠═f3679948-0c5e-4f0a-888b-db6aa614f298
# ╠═1b78fd90-edd0-4225-b6fe-406ba7a4d030
# ╠═859cd515-0b34-460a-bf08-33e43f743695
# ╠═c3650a56-7479-43c0-be7d-4d31d6ce1f06
# ╠═95d7c96c-9a8b-40d1-871a-e3468ad113ce
# ╟─cceffc6a-458f-47fa-832a-6f5711a706b1
# ╠═31091273-f1d4-4a67-935f-2526e0841d98
# ╠═2e379c61-b85f-4ace-bc1a-6b04c6452b56
# ╟─bec867c7-2227-4a5f-aff1-b0a1bc65b398
# ╠═4fdd5a9f-cc78-40a1-bd4b-343e964c4a1f
# ╠═eb42e483-5c95-42e9-9211-a47fedfdadd0
# ╠═0142f9a5-a28f-4f70-8207-f1300ba86c5d
# ╠═6fac34c4-2711-4105-b3df-bf063a166b3b
# ╠═3ed4009b-0f82-470c-b8ec-df3183f7fd66
# ╠═d3102841-af94-4f70-a90d-0f21d80d193e
# ╠═c38b74ce-f51a-4508-9f71-15a11135aefc
# ╠═56438479-98e8-4793-a2f9-6bbc555b7993
# ╠═8c168dc9-6a16-4baf-a2c7-20d7bdc59772
# ╠═15e5185d-a1a7-46e7-a9da-1650e72b3ad0
# ╠═b0bc94b4-f65c-41fb-a72a-bf67cf1c6076
# ╠═016eced6-397f-4c68-a0d8-8de03298dffb
# ╠═d5454c4e-748c-4888-9102-eb6cefbc8409
# ╠═4061e3bf-0b1c-48b5-abf2-76dd3b371f52
# ╠═6980e00b-8bb3-4486-94b7-43095a91c4eb
# ╠═a490b4b1-c21b-41bb-8915-3b2e215df15c
# ╠═81c75b21-6910-47fa-988d-8f240d9e1473
