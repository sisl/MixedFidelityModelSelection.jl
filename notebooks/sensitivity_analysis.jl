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

# ╔═╡ 89f4b7b4-8406-4072-ab3a-a796e0932e42
begin
	using Revise
	using POMDPs
	using POMCPOW
	using Plots; default(fontfamily="Computer Modern", framestyle=:box)
	using Statistics
	using Random
	using Distributions

	using Pkg
	Pkg.develop(path="../../MineralExploration/")
	using MineralExploration
	Pkg.develop(path="../../MixedFidelityModelSelection/")
	using MixedFidelityModelSelection
	
	plot_ore_map = MineralExploration.plot_ore_map
	plot_mass_map = MineralExploration.plot_mass_map
end

# ╔═╡ ebcf80c9-8710-4b2e-b9d0-ed2177a9db5c
using AverageShiftedHistograms

# ╔═╡ f3100c77-7311-4479-b9ca-039f55ff2792
using Images

# ╔═╡ 56a0b3cf-0d61-446b-b83c-d2fd42ef0cd5
using PlutoUI

# ╔═╡ 657bc4a8-fe99-43fe-b258-de271062dc53
md"""
# Sensitivity Analysis
"""

# ╔═╡ 46dbd5fc-12e8-406d-93e4-0bb2e32ca05d
TableOfContents()

# ╔═╡ 0842c42b-b85f-4251-a43b-2b58d958cd91
begin
	N_INITIAL = 0
	MAX_BORES = 25
	MIN_BORES = 5
	GRID_SPACING = 0
	MAX_MOVEMENT = 20
	SAVE_DIR = "./data/sensitivity_analysis/"
	!isdir(SAVE_DIR) && mkdir(SAVE_DIR)
end

# ╔═╡ a6fecc38-adff-4be0-a777-ff85af426dd9
grid_dims = (50, 50, 1)

# ╔═╡ c33957b1-6adc-4d3e-8241-f35d00c09ed1
begin
	true_mainbody = BlobNode(grid_dims=grid_dims)
	mainbody = CircleNode(grid_dims=grid_dims)
end

# ╔═╡ 9d4d4bf0-b03a-4c97-9252-abe85ecd8e4e
# begin
# 	prev_shape_types = ["BlobNode", "EllipseNode", "CircleNode"]
# 	prev_mass_results = generate_ore_mass_samples(grid_dims, prev_shape_types;
# 		N=10_000, apply_scale=false)
# 	# 	true_mainbody = SingleFixedNode(grid_dims=grid_dims)
# 	# 	true_mainbody = MultiVarNode(grid_dims=grid_dims)
# end

# ╔═╡ c36f7e41-2843-4c04-bcf8-a4e760a4a5e2
begin
	shape_types =
		["BlobNode", "EllipseNode", "CircleNode", "SingleFixedNode", "MultiVarNode"]
	# shape_types = [shape_types[2]]
	mass_results, mass_params = generate_ore_mass_samples(grid_dims, shape_types;
		N=1000, apply_scale=true)
end

# ╔═╡ 690e10b0-cea3-4569-a35e-7bb72ddaa0fe
MineralExploration.get_prescaled_parameters(EllipseNode, grid_dims)

# ╔═╡ c70cb301-84b2-48fe-bc8b-a505086886bf
begin
	labels = ["bezier", "ellipse", "circle", #]
		"gaussian (single)", "gaussian (multi)"]
	# labels = [labels[2]]
	dict_mass_results = mass_results
	ore_masses = [
		dict_mass_results["BlobNode"],
		dict_mass_results["EllipseNode"],
		dict_mass_results["CircleNode"],
		dict_mass_results["SingleFixedNode"],
		dict_mass_results["MultiVarNode"],
	]
	plot_ore_mass_distributions(ore_masses, labels, grid_dims)
end

# ╔═╡ d735b9b3-edde-4d74-99a5-5885c9cd0324
begin
	m_ash = 40
	μσ_label = d -> string(round(mean(d), digits=2), " ± ", round(std(d), digits=2))
	Plots.plot(ash(mass_results["BlobNode"], m=m_ash), hist=false,
		label="blob ($(μσ_label(mass_results["BlobNode"])))", c=:crimson)
	Plots.plot!(ash(mass_results["EllipseNode"], m=m_ash), hist=false,
		label="ellipse ($(μσ_label(mass_results["EllipseNode"])))", c=:blue)
	Plots.plot!(ash(mass_results["CircleNode"], m=m_ash), hist=false,
		label="circle ($(μσ_label(mass_results["CircleNode"])))", c=:green)
	title!("ore mass distributions")
end

# ╔═╡ 9ca167dc-fff8-45cd-82f5-6fc42937ebc6
begin
	bins = 20
	histogram(mass_results["BlobNode"], label="blob", alpha=0.3, bins=bins, c=:crimson)
	histogram!(mass_results["EllipseNode"], label="ellipse", alpha=0.3, bins=bins, c=:blue)
	histogram!(mass_results["CircleNode"], label="circle", alpha=0.3, bins=bins, c=:green)
	title!("ore mass distributions")
end

# ╔═╡ 4a861abb-7a67-439e-b166-dd52c265f590
md"""
- Blob: $35.29 \pm 26.25$
- Ellipse: $38.95 \pm 30.25$
- Circle: $49.27 \pm 42.77$
---
- Blob: $36.00 \pm 26.34$
- Ellipse: $66.08 \pm 30.35$
- Circle: $70.64 \pm 39.65$
"""

# ╔═╡ 7520e45e-c5b1-40b4-b7ae-ca95c4566234
md"""
## Standardization
"""

# ╔═╡ 27515004-9203-478d-9e7f-78493f254782
(μᵣ,σᵣ) = 10, 5

# ╔═╡ 59612c9d-8e89-48f8-acd8-c73ec865e171
D = μᵣ .+ σᵣ*randn(10_000)

# ╔═╡ 1d8fdca5-cb07-4740-b39d-5a14a1d2c721
histogram(D)

# ╔═╡ 0396880b-567a-42d8-9740-95f4fe87ea79
function μhistogram(X; kwargs...)
	ρ = x->round(x; digits=2)
	histogram(X, label="$(ρ(mean(X))) ± $(ρ(std(X)))"; kwargs...)
end

# ╔═╡ e936350a-8bc5-41db-88a8-496eb434223d
μhistogram(standardize(D; target_μ=150, target_σ=50))

# ╔═╡ 44d787af-0dc7-40a2-a75a-a6147c71810f
μ, σ = calculate_standardization(D)

# ╔═╡ f277cce5-a9c8-4a3e-b921-9797ba6a6fe9
μhistogram(standardize(D, μ, σ; target_μ=150, target_σ=50))

# ╔═╡ fa42dadc-c246-469b-a847-e79619728395
@bind x PlutoUI.Slider(0:20, default=10, show_value=true)

# ╔═╡ 0a5b69bd-de40-494e-929a-b9429c1cb420
standardize(x, μ, σ; target_μ=150, target_σ=50)

# ╔═╡ 9523ea5c-1e82-4ff3-b9a9-c424dc467b2c
md"""
# Truth observation (high-fidelity)
"""

# ╔═╡ 4f701ab7-a8da-4465-8d09-ca92bcc39cf6
function high_fidelity_obs(m, subsurface_map, a)
	hf_coords = trunc.(Int, a.coords ./ m.ratio[1:2])
	return subsurface_map[hf_coords[1], hf_coords[2], 1]
end

# ╔═╡ 2ea1631f-47ef-49c0-92c7-9ee9d09501d8
@bind high_fidelity CheckBox(true)

# ╔═╡ dc115277-7458-413a-9e2e-a166cbca069e


# ╔═╡ 4cbaf630-1793-4250-a3ef-e6a59d33d8df
md"""
# Normal samples of ore mass
"""

# ╔═╡ 56b34fe6-91f0-45fe-9d62-aeaff821f619
md"""
---
"""

# ╔═╡ 836c9078-5899-439b-b62c-9bb14fe187ea
xydim = 10

# ╔═╡ 7597e5aa-3fce-484f-8667-1356d086b6a5
begin
	m = MineralExplorationPOMDP(max_bores=MAX_BORES, delta=GRID_SPACING+1,
		                        grid_spacing=GRID_SPACING,
	                            true_mainbody_gen=true_mainbody,
		                        mainbody_gen=mainbody,
		                        original_max_movement=MAX_MOVEMENT,
	                            min_bores=MIN_BORES, grid_dim=(xydim,xydim,1))
	initialize_data!(m, N_INITIAL)
end

# ╔═╡ e71d1a1f-3f81-46a5-b2cf-a95d4c9f4d57
D_masses, D_mass_samples, _ = generate_ore_mass_samples(m);

# ╔═╡ ab663118-b12f-4bbe-bcba-a90491af4056
begin
	i_mass = 1
	D_ore = D_mass_samples[i_mass][:,:,1]
end;

# ╔═╡ b0366fee-af01-4a61-be8a-7be97cae9c76
begin
	super_ore = D_mass_samples[argmax(D_masses)]
	Plots.plot(plot_ore_map(super_ore),
		   first(plot_mass_map(m, super_ore)),
		   size=(700,300))
end

# ╔═╡ 9b2abacb-2ce5-4590-978d-a8827b12091c
begin
	Random.seed!(2) # Blob 12 (2 10x10), Ellipse 12, Circle 2
	masses, ore_samples, ore_params = generate_ore_mass_samples(m;
		N=4, apply_scale=true)
end;

# ╔═╡ 9239edb5-7607-4d12-9d93-edd34c0094d8
begin

	ore_sample_figs = []
	for i in 1:length(ore_samples)
		ore_sample = ore_samples[i]
		push!(ore_sample_figs, plot_ore_map(ore_sample))
		push!(ore_sample_figs, first(plot_mass_map(m, ore_sample)))
	end
	Plots.plot(ore_sample_figs...; cbar=false,
			   size=(700,300*length(ore_samples)),
		       layout=@layout([a b; c d; e f; g h]))
end

# ╔═╡ 81c94714-bdb4-4005-81f8-4332acea48b2
begin
	seed = 1994358009 # 1268077864
	Random.seed!(seed) # determinism
	ds0 = POMDPs.initialstate_distribution(m)
	s0 = rand(ds0; truth=true)
end;

# ╔═╡ 3dea8e40-bb95-4ea6-8db9-4576fdf09a1c
s0.mainbody_map |> plot_ore_map

# ╔═╡ 3c3eddc9-bfb7-4af8-98ba-f583a95126ec
s0.ore_map |> plot_ore_map; title!("true (noisy) ore field")

# ╔═╡ b2cebfa0-59a6-4af5-9997-610ab1e5db90
begin
	Random.seed!(0)
	s = rand(ds0; truth=true)
	lower_dims = (10,10,1)
	s_mainbody_map_scaled = imresize(s.mainbody_map[:,:,1], lower_dims[1:2])
	m.dim_scale
end

# ╔═╡ 71f6c5b2-e21d-42dd-928f-f970746267e2
s.ore_map |> plot_ore_map

# ╔═╡ daa61af3-8b65-4822-9551-afb2ef01182f
ratio = lower_dims ./ (50, 50, 1)

# ╔═╡ 9a139e02-a6f6-4158-97cb-9a7e32546cae
dim_scale = 1/prod(ratio)

# ╔═╡ d061092d-504b-4f88-956e-9aae058b9c38
@bind xᵢ Slider(1:lower_dims[1], default=lower_dims[1]÷2)

# ╔═╡ 7e03b9d1-de05-45c6-a5a4-606a58a4bc67
@bind yᵢ Slider(1:lower_dims[2], default=lower_dims[2]÷2)

# ╔═╡ bc54910a-12a3-4670-96fe-24de45084bae
a = [xᵢ,yᵢ]

# ╔═╡ 4216d0b6-65eb-4257-adfc-2607b7d2f82d
a ./ ratio[1:2]

# ╔═╡ 265f921d-973d-4d11-9b29-20c56094df5d
index = trunc.(Int, a ./ ratio[1:2]) # round?

# ╔═╡ 2c0a4494-5ba0-40f7-bf62-d74d5816dcd0
index_round = round.(Int, a ./ ratio[1:2]) # round?

# ╔═╡ d129bf95-5e51-4f9d-9de4-4c234084973d
A = (coords=a,)

# ╔═╡ 0148a0d2-e5bf-43a6-bfe0-26533372c44d
obs_high_fidelity_trunc = s.mainbody_map[index[1], index[2], 1]

# ╔═╡ 9b7216b2-7b4d-44ba-acf6-1cbe8907cfd0
obs_high_fidelity_round = s.mainbody_map[index_round[1], index_round[2], 1]

# ╔═╡ 23a060b5-f0fe-4a9c-b4e3-6be5f074e60b
obs_low_fidelity = s_mainbody_map_scaled[a[1], a[2]]

# ╔═╡ b472565b-46a5-4bea-8562-f964f38c6547
abs(obs_high_fidelity_trunc - obs_low_fidelity)

# ╔═╡ 486f29e1-8905-4fd7-b37a-6ab727479e51
abs(obs_high_fidelity_round - obs_low_fidelity)

# ╔═╡ 200dc66a-fdcc-46f7-8123-ef4153431831
M = (ratio=lower_dims ./ m.grid_dim, )

# ╔═╡ 400e3af1-765b-4984-9d1c-6324b1a77098
high_fidelity_obs(M, s.mainbody_map, A)

# ╔═╡ a077bacf-621c-43de-b9d6-a2ba7eb975da
high_fidelity_obs(M, s.ore_map, A)

# ╔═╡ 4aea7b45-2486-4249-acbb-202d3af78afb
begin
	err_trunc = []
	err_round = []
	for xᵢ in 1:lower_dims[1]
		for yᵢ in 1:lower_dims[2]
			a = [xᵢ,yᵢ]
			index_trunc = trunc.(Int, a ./ ratio[1:2])
			index_round = round.(Int, a ./ ratio[1:2])
			obs_high_fidelity_trunc = s.mainbody_map[index_trunc[1],index_trunc[2],1]
			obs_high_fidelity_round = s.mainbody_map[index_round[1],index_round[2],1]
			obs_low_fidelity = s_mainbody_map_scaled[a[1], a[2]]
			push!(err_trunc, abs(obs_high_fidelity_trunc - obs_low_fidelity))
			push!(err_round, abs(obs_high_fidelity_round - obs_low_fidelity))
		end
	end
	err_trunc, err_round
end

# ╔═╡ ca04f24f-921f-4e77-9a79-8cf73a406d76
mean(err_trunc), std(err_trunc)

# ╔═╡ 3eaf3478-e970-4d58-8243-911b0a7ef65b
mean(err_round), std(err_round)

# ╔═╡ 602affac-25ef-4f0c-aa8c-d34958bc4ff7
if high_fidelity
	s.mainbody_map |> plot_ore_map
	scatter!([index[2]], [index[1]], marker=:square, c=:red, label="obs")
else
	s_mainbody_map_scaled |> plot_ore_map
	scatter!([a[2]], [a[1]], c=:red, marker=:square, label="obs")
end

# ╔═╡ e45a4fe0-b82d-4101-aa77-eff2a4aaf941
gp_ore_map = Base.rand(ds0.rng, ds0.gp_distribution, 1); plot_ore_map(gp_ore_map)

# ╔═╡ 1ab44e34-b0d4-477a-90a0-9f578e1eedd9
begin
	plot_ore_map(imresize(s0.ore_map[:,:,1], (xydim,xydim)))
	title!("true ore field $(xydim)x$(xydim)")
end	

# ╔═╡ 26638f04-b7d2-49a1-9d20-a774c7ef232e
begin
	mass_fig, r_massive = plot_mass_map(m, s0.ore_map)
	mass_fig
end

# ╔═╡ d2dfed3e-108d-4c92-b1ae-c0acc04d3b87
#   # Check the difference!
# 	@info (r_massive_before, r_massive, (r_massive - A*r_massive_before))

# ╔═╡ a53e6f08-5886-4f14-a24a-95d02837ba09
md"""
## Sampling massive ore bodies (seeds)
"""

# ╔═╡ eb4006dd-3070-4ea0-951c-f5d2dec7874f
begin
	Random.seed!(1)
	massive_samples, massive_seeds = sample_massive(m, [150, 160]; N=1, max_iters=1000)
end

# ╔═╡ f0251439-ee6c-4fe5-b5af-624364d2e2fd
massive_seeds

# ╔═╡ 54355228-09b9-40ee-b108-66dc261119da
Plots.plot([plot_mass_map(m, sample)[1] for sample in massive_samples]...)

# ╔═╡ 8c99bca8-aa9f-4e8d-8857-bf3530d8160e
Plots.plot([plot_ore_map(sample) for sample in massive_samples]...)

# ╔═╡ 6e67449f-7690-4dc9-bf94-b26e01a4c033
md"""
# Planner
"""

# ╔═╡ 0da45e59-e343-40f1-bb00-1bf0feaf559e
begin
	up = MEBeliefUpdater(m, 1000, 2.0)
	b0 = POMDPs.initialize_belief(up, ds0)
end;

# ╔═╡ 480cf970-0dfc-4aae-acd7-993df09eecfd
begin
	next_action = NextActionSampler()
	tree_queries = [100, 1_000, 10_000]
	i_tree_queries = 1
	solver = POMCPOWSolver(tree_queries=tree_queries[i_tree_queries],
	                       check_repeat_obs=true,
	                       check_repeat_act=true,
	                       next_action=next_action,
	                       k_action=2.0,
	                       alpha_action=0.25,
	                       k_observation=2.0,
	                       alpha_observation=0.1,
	                       criterion=POMCPOW.MaxUCB(100.0),
	                       final_criterion=POMCPOW.MaxQ(),
	                       # final_criterion=POMCPOW.MaxTries(),
	                       estimate_value=0.0
	                       # estimate_value=leaf_estimation
	                       )
	planner = POMDPs.solve(solver, m)
end

# ╔═╡ f17a096e-615b-4525-967f-acff481bcd2f
@bind perform_run CheckBox()

# ╔═╡ 23adfaf7-2bf9-43fd-b699-78c12615fda0
trial_seed = 0

# ╔═╡ 0acbcbd8-54ef-4292-af97-9263c669dab2
function run_single_trial(seed)
	Random.seed!(seed)
	timing = @timed results =
		run_trial(m, up, planner, s0, b0, save_dir=nothing, display_figs=false)
	return results, timing
end

# ╔═╡ 3f6727e4-9edb-4565-b9f7-2c5b4711c38b
# if perform_run
# 	results, timing = run_single_trial(trial_seed)
# end

# ╔═╡ 25fa502b-6655-47e4-a659-587abcba5129
# if perform_run
# 	timing.time
# end

# ╔═╡ 6985accd-9080-4a90-a92d-a704525fc5b3
# if perform_run
# 	abs_errs, rel_errs = results[3:4]
# end

# ╔═╡ e294f086-6169-48f6-a40f-cda9c4907d86
# if perform_run
#     ts = [1:length(rel_errs);] .- 1
#     plot(ts, rel_errs, title="relative volume error",
#         xlabel="time step", ylabel="relative error", legend=:none, lw=2, c=:crimson)
# 	rel_err_fig = plot!([xlims()...], [0, 0], lw=1, c=:gray, xlims=xlims())
# end

# ╔═╡ 019bae5d-37fb-4d2b-b8a3-9ec6531ca9b7
r1 = first.(poutput)

# ╔═╡ b2039618-72b4-415a-a280-e5a2d6eba5da
r2 = last.(poutput)

# ╔═╡ Cell order:
# ╟─657bc4a8-fe99-43fe-b258-de271062dc53
# ╠═89f4b7b4-8406-4072-ab3a-a796e0932e42
# ╠═46dbd5fc-12e8-406d-93e4-0bb2e32ca05d
# ╠═0842c42b-b85f-4251-a43b-2b58d958cd91
# ╠═a6fecc38-adff-4be0-a777-ff85af426dd9
# ╠═c33957b1-6adc-4d3e-8241-f35d00c09ed1
# ╠═9d4d4bf0-b03a-4c97-9252-abe85ecd8e4e
# ╠═c36f7e41-2843-4c04-bcf8-a4e760a4a5e2
# ╠═690e10b0-cea3-4569-a35e-7bb72ddaa0fe
# ╠═c70cb301-84b2-48fe-bc8b-a505086886bf
# ╠═ebcf80c9-8710-4b2e-b9d0-ed2177a9db5c
# ╠═d735b9b3-edde-4d74-99a5-5885c9cd0324
# ╠═9ca167dc-fff8-45cd-82f5-6fc42937ebc6
# ╟─4a861abb-7a67-439e-b166-dd52c265f590
# ╟─7520e45e-c5b1-40b4-b7ae-ca95c4566234
# ╠═27515004-9203-478d-9e7f-78493f254782
# ╠═59612c9d-8e89-48f8-acd8-c73ec865e171
# ╠═1d8fdca5-cb07-4740-b39d-5a14a1d2c721
# ╠═0396880b-567a-42d8-9740-95f4fe87ea79
# ╠═e936350a-8bc5-41db-88a8-496eb434223d
# ╠═44d787af-0dc7-40a2-a75a-a6147c71810f
# ╠═f277cce5-a9c8-4a3e-b921-9797ba6a6fe9
# ╠═fa42dadc-c246-469b-a847-e79619728395
# ╠═0a5b69bd-de40-494e-929a-b9429c1cb420
# ╟─9523ea5c-1e82-4ff3-b9a9-c424dc467b2c
# ╠═e71d1a1f-3f81-46a5-b2cf-a95d4c9f4d57
# ╠═b0366fee-af01-4a61-be8a-7be97cae9c76
# ╠═ab663118-b12f-4bbe-bcba-a90491af4056
# ╠═3dea8e40-bb95-4ea6-8db9-4576fdf09a1c
# ╠═3c3eddc9-bfb7-4af8-98ba-f583a95126ec
# ╠═71f6c5b2-e21d-42dd-928f-f970746267e2
# ╠═f3100c77-7311-4479-b9ca-039f55ff2792
# ╠═b2cebfa0-59a6-4af5-9997-610ab1e5db90
# ╠═daa61af3-8b65-4822-9551-afb2ef01182f
# ╠═9a139e02-a6f6-4158-97cb-9a7e32546cae
# ╠═bc54910a-12a3-4670-96fe-24de45084bae
# ╠═4216d0b6-65eb-4257-adfc-2607b7d2f82d
# ╠═d061092d-504b-4f88-956e-9aae058b9c38
# ╠═7e03b9d1-de05-45c6-a5a4-606a58a4bc67
# ╠═265f921d-973d-4d11-9b29-20c56094df5d
# ╠═2c0a4494-5ba0-40f7-bf62-d74d5816dcd0
# ╠═0148a0d2-e5bf-43a6-bfe0-26533372c44d
# ╠═9b7216b2-7b4d-44ba-acf6-1cbe8907cfd0
# ╠═23a060b5-f0fe-4a9c-b4e3-6be5f074e60b
# ╠═b472565b-46a5-4bea-8562-f964f38c6547
# ╠═486f29e1-8905-4fd7-b37a-6ab727479e51
# ╠═4f701ab7-a8da-4465-8d09-ca92bcc39cf6
# ╠═400e3af1-765b-4984-9d1c-6324b1a77098
# ╠═a077bacf-621c-43de-b9d6-a2ba7eb975da
# ╠═200dc66a-fdcc-46f7-8123-ef4153431831
# ╠═d129bf95-5e51-4f9d-9de4-4c234084973d
# ╠═4aea7b45-2486-4249-acbb-202d3af78afb
# ╠═ca04f24f-921f-4e77-9a79-8cf73a406d76
# ╠═3eaf3478-e970-4d58-8243-911b0a7ef65b
# ╠═2ea1631f-47ef-49c0-92c7-9ee9d09501d8
# ╠═602affac-25ef-4f0c-aa8c-d34958bc4ff7
# ╠═dc115277-7458-413a-9e2e-a166cbca069e
# ╟─4cbaf630-1793-4250-a3ef-e6a59d33d8df
# ╠═e45a4fe0-b82d-4101-aa77-eff2a4aaf941
# ╠═9b2abacb-2ce5-4590-978d-a8827b12091c
# ╠═9239edb5-7607-4d12-9d93-edd34c0094d8
# ╟─56b34fe6-91f0-45fe-9d62-aeaff821f619
# ╠═836c9078-5899-439b-b62c-9bb14fe187ea
# ╠═7597e5aa-3fce-484f-8667-1356d086b6a5
# ╠═81c94714-bdb4-4005-81f8-4332acea48b2
# ╠═1ab44e34-b0d4-477a-90a0-9f578e1eedd9
# ╠═26638f04-b7d2-49a1-9d20-a774c7ef232e
# ╠═d2dfed3e-108d-4c92-b1ae-c0acc04d3b87
# ╟─a53e6f08-5886-4f14-a24a-95d02837ba09
# ╠═eb4006dd-3070-4ea0-951c-f5d2dec7874f
# ╠═f0251439-ee6c-4fe5-b5af-624364d2e2fd
# ╠═54355228-09b9-40ee-b108-66dc261119da
# ╠═8c99bca8-aa9f-4e8d-8857-bf3530d8160e
# ╟─6e67449f-7690-4dc9-bf94-b26e01a4c033
# ╠═0da45e59-e343-40f1-bb00-1bf0feaf559e
# ╠═480cf970-0dfc-4aae-acd7-993df09eecfd
# ╠═56a0b3cf-0d61-446b-b83c-d2fd42ef0cd5
# ╠═f17a096e-615b-4525-967f-acff481bcd2f
# ╠═23adfaf7-2bf9-43fd-b699-78c12615fda0
# ╠═0acbcbd8-54ef-4292-af97-9263c669dab2
# ╠═3f6727e4-9edb-4565-b9f7-2c5b4711c38b
# ╠═25fa502b-6655-47e4-a659-587abcba5129
# ╠═6985accd-9080-4a90-a92d-a704525fc5b3
# ╠═e294f086-6169-48f6-a40f-cda9c4907d86
# ╠═019bae5d-37fb-4d2b-b8a3-9ec6531ca9b7
# ╠═b2039618-72b4-415a-a280-e5a2d6eba5da
