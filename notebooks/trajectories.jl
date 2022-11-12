### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ a333c86e-d90e-49e5-8623-8e844ad2a9b8
begin
	using Revise
	using Pkg
	Pkg.develop(path="..//..//MineralExploration//")
	using MineralExploration
	Pkg.develop(path="..//")
	using MixedFidelityModelSelection
	Pkg.develop(path="..//scripts//MEParallel.jl//")
	using MEParallel
	using BSON
	using Statistics
	using PlutoUI
	using Plots; default(fontfamily="Computer Modern", framestyle=:box)
end

# ╔═╡ 74555820-8c6e-4714-86f2-d6c33434bb80
using Reel

# ╔═╡ deddcc30-62e4-4736-830b-550d880970d1
using Distributions

# ╔═╡ 4bec3cad-4e8f-4b59-b878-b584cf8daacd
using Random, LinearAlgebra

# ╔═╡ 1c6b50e0-7024-45a2-aa6b-02031eb880a4
using PDMatsExtras

# ╔═╡ fb8868be-997f-11ec-0fb4-37f4c70c68c7
md"""
# Trajectories
"""

# ╔═╡ 62013d0f-8cbd-473a-8ed8-9662a567328b
results_regret = BSON.load("..\\scripts\\MEParallel.jl\\results\\results_500seeds.bson")[:results]

# ╔═╡ 7d875d7d-e2e8-413f-aece-452b539fca16
shapekeys = [:blob, :ellipse, :circle]

# ╔═╡ ba7a7f30-fc84-4a5e-b853-9946c69a1920
function fargmin(results, f)
	keyset = collect(keys(results))
	return keyset[argmin(f(results[k]) for k in keyset)]
end

# ╔═╡ d1c80a9e-e5a9-458a-b76f-8dfaa5ab8b53
function fargmax(results, f)
	keyset = collect(keys(results))
	return keyset[argmax(f(results[k]) for k in keyset)]
end

# ╔═╡ fe64f46b-98f6-454a-9ab8-2f73a9722e11
mean_regret = res->mean(regret(res))

# ╔═╡ 3ac672e1-33de-4405-ace2-b982ab8f94df
quant_regret = α->(res->quantile(regret(res), α)) # Currying

# ╔═╡ 4185c93a-42b1-43b8-b44d-7c301b270392
md"""
###### Minimum accuracy
"""

# ╔═╡ 77ca5cd5-1501-48d9-a9ba-f3b0c482c8f2
τ_min_accuracy = (:circle, (30, 30, 1), 100) # fargmin(results_regret, accuracy)

# ╔═╡ 06a55c95-92cc-453f-8e10-1d53f153fb23
md"""
###### Maximum expected regret
"""

# ╔═╡ 4a144a76-625c-4b17-b998-89b7fef5014d
τ_max_expected_regret = fargmax(results_regret, mean_regret)

# ╔═╡ 80bbd217-d81e-4b80-b455-a84309c4fe8d
md"""
###### Maximum 90th percentile of regret
"""

# ╔═╡ fc372d8f-2f97-4f5c-8062-c35f7864912f
τ_max_90th_regret = fargmax(results_regret, quant_regret(0.9))

# ╔═╡ fa227921-0730-4f33-ae61-8ea635be9e44
md"""
###### Maximum bias
"""

# ╔═╡ 4b00b8cd-b58c-40a1-b86b-adf77918ec3c
f(res, extraction_cost=150) = res[:r_massive]

# ╔═╡ ff7463c2-1f8f-44ae-85e0-fff247bd5560
f̂(res) = f(res) .+ last.(res[:rel_errors])

# ╔═╡ ed628d8b-cf6e-4e4a-865b-d844f676d512
bias(res) = mean(mean(f̂(res)) .- f(res))

# ╔═╡ 916b42b0-b171-44fc-b27c-960208c2f843
τ_max_bias = fargmax(results_regret, bias)

# ╔═╡ 7a23a34c-a0e4-47dd-970f-f42474bde9c5
md"""
# Bias
"""

# ╔═╡ 321013a9-6371-4628-adaf-ca89017041e0
wc_bias_res = results_regret[τ_max_bias]

# ╔═╡ 9f32913d-7913-47f5-90f6-dbfa6955e46d
bias(wc_bias_res)

# ╔═╡ 1ce9a2c4-917e-47e7-bb0b-4000aed7ae3d
wc_bias_idx = 7 # 8, 13

# ╔═╡ dd1fab4d-1c10-4a55-942b-45cdc3a3b96f
wc_bias_res[:r_massive]

# ╔═╡ 34d07754-d716-4c4d-a346-56b73947b965
last.(wc_bias_res[:rel_errors])[wc_bias_idx]

# ╔═╡ 56f2a8ec-192a-4d5b-ae82-376062bbac8e
wc_bias_res[:r_massive][wc_bias_idx]

# ╔═╡ d4a91997-9e5d-4b3e-809d-f0532f06e07d
wc_bias_config = results_regret[τ_max_bias][:config][wc_bias_idx]

# ╔═╡ 3a33df07-bcff-4f22-abb8-565613cd8fa9
# results_bias, final_belief, initial_belief = evaluate(trial_bias; save_dir="worst_case_bias")

# ╔═╡ 729baffb-154b-4e30-8918-8547c1071177
function create_mass_gif(particles; output="particles_mass.gif", fps=3)
	frames = Frames(MIME("image/png"), fps=fps)

	local frame
	for particle in particles
		frame = MineralExploration.plot_mass_map(particle.ore_map,0.7)[1]
		push!(frames, frame)
	end
	write(output, frames)
	LocalResource("./$output")
end

# ╔═╡ 9c1c3ebf-a63a-41c1-9f7d-6d853e282a74
# create_mass_gif(initial_belief.particles)

# ╔═╡ 776b4f28-81ff-44ef-b6c2-8b418075cc1e
initial_belief.particles[1].ore_map |> x->MineralExploration.plot_mass_map(x,0.7)[1]

# ╔═╡ 88c4fc71-fb7d-4b44-9cfd-77c100437c3a
MineralExploration.summarize(initial_belief)[1] |> x->MineralExploration.plot_mass_map(x,0.7)[1]

# ╔═╡ 1d4e6d83-e96b-4138-b6c1-844b20fe3729
begin
	MineralExploration.summarize(initial_belief)[1] |> MineralExploration.plot_ore_map
	title!("initial belief")
end

# ╔═╡ e74da8ce-ff04-4778-9ec6-25ce30ed7cd2
MineralExploration.summarize(final_belief)[1] |> x->MineralExploration.plot_mass_map(x,0.7)[1]

# ╔═╡ f253b3b3-bcac-4330-af1a-6cd098d4375f
MineralExploration.summarize(final_belief)[1] |> MineralExploration.plot_ore_map

# ╔═╡ 6b53029d-ffe8-4fb8-ac97-9b1b5677821e
results_bias_dict = MEParallel.convert(Dict, results_bias);

# ╔═╡ 7bc1f9c8-b771-4538-8970-05c4d4f5d12f
results_bias_dict

# ╔═╡ 44a5c6ca-59de-4577-82c5-8fb87cdb6a47
single_bias(res) = res[:rel_errors][end]

# ╔═╡ 4f6d7b0a-e424-4348-a566-5300a0ce1420
single_bias(results_bias_dict)

# ╔═╡ 8daff9ab-c4fc-4d34-95fc-a47ab70aec50
max_biases = last.(results_regret[τ_max_bias][:rel_errors])

# ╔═╡ b15ad0b4-37c9-4533-9d1c-5bc797079043
# volumes = last

# ╔═╡ 71782ce8-1fdf-4140-a53a-9ad70ce5779f
md"""
# Worst-case: Accuaracy
"""

# ╔═╡ 506dc301-81e8-4c80-bf7e-3843f9b911d2
τ_min_accuracy

# ╔═╡ 6acd1e62-ca92-4f62-96a6-38bc71ec6871
wc_acc_res = results_regret[τ_min_accuracy]

# ╔═╡ 2f03b56f-4ea7-4dd3-8974-ecd570ff6c0b
wc_acc_idx = 18 # 18 # argmax(regret(results_regret[τ_min_accuracy]))

# ╔═╡ 6b756e10-0f71-4ccb-bc23-ff9a74bba135
regret(results_regret[τ_min_accuracy])[wc_acc_idx]

# ╔═╡ 2a7d356e-20f8-403b-9e9b-726e54919643
wc_acc_config = results_regret[τ_min_accuracy][:config][wc_acc_idx]

# ╔═╡ b850b7fc-641e-498d-87ff-43f6e938abb6
function convert(::Type{MEConfiguration}, config::Dict)
	seed = config[:seed]
	grid_dims = config[:grid_dims]
	pomcpow_iters = config[:pomcpow_iters]
	mainbody_type = eval(Meta.parse(config[:mainbody_type]))
	params = convert(MEJobParameters, config[:params])
	return MEConfiguration(seed, grid_dims, pomcpow_iters, mainbody_type, params)
end

# ╔═╡ d28a33db-b17c-415d-923b-8cdf87f5c860
function convert(::Type{MEJobParameters}, params::Dict)
	high_fidelity_dim = params[:high_fidelity_dim]
	min_bores = params[:min_bores]
	max_bores = params[:max_bores]
	grid_spacing = params[:grid_spacing]
	max_movement = params[:max_movement]
	n_initial = params[:n_initial]
	name = params[:name]
	return MEJobParameters(high_fidelity_dim=high_fidelity_dim,
                           min_bores=min_bores,
                           max_bores=max_bores,
                           grid_spacing=grid_spacing,
                           max_movement=max_movement,
                           n_initial=n_initial,
                           name=name)
end

# ╔═╡ ebb31bbf-cee0-4a09-a044-b864e0856e84
config_bias = convert(MEConfiguration, wc_bias_config)

# ╔═╡ 816a3b19-8e43-4b84-86fc-893c4ed92619
trial_bias = initialize(config_bias)

# ╔═╡ de18f869-8f5e-4c75-b631-da8ae0f61794
config = convert(MEConfiguration, wc_acc_config)

# ╔═╡ eb9cd498-27b6-4747-9d84-5b02b73c9eae
trial = initialize(config)

# ╔═╡ c781d26c-bbcc-4ea9-a2d4-128d3dfd9bd4
results = evaluate(trial; save_dir="worst_case_accuracy")

# ╔═╡ 235d6716-2908-497f-9001-4b1ed65b5a81
results_dict = MEParallel.convert(Dict, results);

# ╔═╡ 9d5c8911-3934-465f-9337-b12ef66ca7b2
regret(results_dict)

# ╔═╡ 6387dd41-dfc4-4109-8602-077dfb13757e
results.ore_map |> MineralExploration.plot_ore_map

# ╔═╡ 2a4125ad-2e2a-4e12-ac38-04ed56c45c2c
md"""
# Worst-case: 90th percentile regret
"""

# ╔═╡ 081c6ba0-4155-4d96-a8d8-b387f2cc5244
wc_regret_res = results_regret[τ_max_90th_regret]

# ╔═╡ a45ba4d6-df96-47c5-a129-a6eab357bed3
wc_regret_res[:config]

# ╔═╡ ae7606ef-fc61-403d-8136-9e1ce23a6a84
md"""
# Solution ideas
"""

# ╔═╡ 28840030-6459-4c6f-a38e-470eb5d9c1c6
fargmin(results_regret, quant_regret(0.9))

# ╔═╡ fbf32114-5237-4997-96c4-cc5475862bdd
fargmin(results_regret, quant_regret(0.95))

# ╔═╡ 56b40041-50fa-43f3-b260-f3b189bf21a7
fargmin(results_regret, mean_regret)

# ╔═╡ 211cbb7d-d2d8-4e5b-9fbe-6a90bc543985
function regret_distribution(results, k)
	histogram(regret(results[k]), bins=5, normalize=:pdf, c=:gray)
	𝒩 = fit(Normal, regret(results[k]))
	𝒩 = TruncatedNormal(𝒩.μ, 𝒩.σ, 0, Inf)
	plot!(x->pdf(𝒩,x), c=:red, lw=2)
	title!(string(k))
end

# ╔═╡ 4c1554ea-7451-4edd-b2d6-d4c56d17d7f0
regret_distribution(results_regret, τ_max_expected_regret)

# ╔═╡ 675f654e-ad55-432a-a421-ca9a28737353
md"""
# Plots
"""

# ╔═╡ e11f7f67-ac9b-4df5-a678-8172b499152f
begin
	MEParallel.Random.seed!(13)
	rand(MEParallel.POMDPs.initialstate_distribution(trial.problem); truth=true).ore_map |> 	MineralExploration.plot_ore_map
	title!("")
	# scatter!([20],[20], label=false, c=:red)
end

# ╔═╡ c96dd779-621a-4c16-b4eb-375b4d67c299
md"""
# Distance to center
"""

# ╔═╡ 763e4c58-72a6-4286-b9ea-c825309105dd
function massive(ore_map, massive_threshold)
	return ore_map .>= massive_threshold
end

# ╔═╡ b8d17852-8668-42ca-98fa-e7f1424ffeb2
function massive_center(mass)
	return mean(findall(mass))
end

# ╔═╡ 55433d68-5dd6-4d1e-a0c7-abde02d7cf65
function plot_massive_center!(center; color=:red, label="center")
	return scatter!([center[2]], [center[1]], c=color, label=label)
end

# ╔═╡ 61106a23-8901-450f-a58a-05c58163f3a8
# used in: mean(Array{CartesianIndex{3}})
function Base.:(/)(c::CartesianIndex{3}, n::Int64)
	return c.I ./ n
end

# ╔═╡ 296c6deb-1ce4-4c9e-af85-2977ffa29990
sum(results_regret[τ_max_bias][:last_action] .== :mine) / 500

# ╔═╡ b8450d7d-ab0f-4429-9a82-169d4eda5fd7
b_mean, _ = MineralExploration.summarize(trial.b0)

# ╔═╡ db0123bc-103e-4ab5-a60f-02bef918dfd9
begin
	MineralExploration.plot_ore_map(b_mean)
	b_max = argmax(b_mean).I
	scatter!([b_max[2]], [b_max[1]], c=:white, label="belief center")
end

# ╔═╡ 19ecdfb4-4cfa-445b-af60-92d376cd7fd1
b_max

# ╔═╡ 6cd5b0cd-685c-4c89-af98-59664a64c95b
b_mean[b_max...]

# ╔═╡ b8a6cf63-04ee-480f-82a4-f2f8a9cb20b1
maximum(b_mean)

# ╔═╡ 8117557b-035d-49bc-801a-11cff86208e3
mean(b_mean)

# ╔═╡ 3ba758ab-d8da-4f67-9b46-b941fe35f4fe
μ = vec(mean(b_mean[:,:,1], dims=2));

# ╔═╡ 2eb69b52-bc1b-469d-8829-de6108d95a01
σ = PSDMat(cov(b_mean[:,:,1], dims=2));

# ╔═╡ da2e51af-1b80-4927-ba7c-b361891662f7
mv = MvNormal(μ, σ)

# ╔═╡ ed2e6025-aa46-4565-8703-e5ad4a8656ce
function mean_distance_to_center(trial::Trial)
	m_thresh = trial.problem.massive_threshold
	true_mass_map = massive(trial.s0.ore_map, m_thresh)
	ratio = trial.problem.ratio
	return mean_distance_to_center(trial.b0, true_mass_map, ratio)
end	

# ╔═╡ 3ec60a22-40b4-4750-88cb-509621aa1c52
function mean_distance_to_center(belief_center, true_mass_map, ratio)
	true_mass_center = massive_center(true_mass_map)
	belief_center_adjusted = belief_center ./ ratio
	return norm(belief_center_adjusted .- true_mass_center)
end

# ╔═╡ 303198db-a517-483f-abc1-ff27f921636a
function belief_center(b0)
	b_mean, _ = MineralExploration.summarize(b0)
	return argmax(b_mean).I
end

# ╔═╡ d32a110e-ca3e-4ba2-8d92-baf93883f67a
function mean_distance_to_center(b0::MineralExploration.MEBelief, true_mass_map, ratio)
	b_center = belief_center(b0)
	return mean_distance_to_center(b_center, true_mass_map, ratio)
end	

# ╔═╡ b0258367-59ce-4b98-b506-f2de7a64f333
mean_distance_to_center(trial)

# ╔═╡ 41ae1023-32c1-4dac-b134-2b01ddf2863d
ratio = (30,30,1) ./ (50,50,1)

# ╔═╡ fe849f40-ace8-41f7-9f10-099209342d67
begin
	true_mass = massive(trial.s0.ore_map, trial.problem.massive_threshold)
	MineralExploration.plot_ore_map(true_mass)
	true_mass_center = massive_center(true_mass)
	plot_massive_center!(true_mass_center)
	plot_massive_center!(b_max ./ ratio;
						 color=:white, label="belief center")
end

# ╔═╡ 2576c794-419b-4de0-8b71-6f5c50a2b9cb
begin
	mass = massive(trial.b0.particles[3].ore_map, trial.problem.massive_threshold)
	MineralExploration.plot_ore_map(mass)
	mass_center = massive_center(mass)
	plot_massive_center!(mass_center)
	plot_massive_center!(true_mass_center .* trial.problem.ratio;
						 color=:white, label="true center")
end

# ╔═╡ 80f8f9ca-1edb-4305-9aac-178bd494c078
mean_distance_to_center(trial.b0, true_mass, (30,30,1) ./ (50,50,1))

# ╔═╡ 1bf3c03a-34e7-4a8e-934a-691afa0586d6
norm((12,12,1) .- true_mass_center .* (30,30,1) ./ (50,50,1))

# ╔═╡ f474bc29-4161-4e66-8ab5-a41f3cdaf670
massive_center(true_mass)

# ╔═╡ 7d41f3be-20d8-4726-b0f9-fd732e4e4abe
b_max ./ ratio

# ╔═╡ f78a4a1a-8075-4393-8736-5f5dfd61c237
norm(massive_center(true_mass) .- b_max ./ ratio)

# ╔═╡ 729e1916-92c6-447c-bfd0-8094820feff8
mean_distance_to_center(b_max, true_mass, ratio)

# ╔═╡ f77a8482-e4d6-47ab-b186-98a2e3a7ed03
b_max ./ ratio

# ╔═╡ 8a1d7ee7-95db-49bd-a9e2-51f17ea07a6a
true_mass_center

# ╔═╡ c2fa1f0f-4ba3-4269-9ff4-be84efc06424
true_mass_center .* trial.problem.ratio

# ╔═╡ ff9561f5-ead3-4849-a97e-805ba844d60a
begin
	b0_set = Dict()
	for shape in [:circle, :ellipse, :blob]
		for grid_dim in [10, 30, 50]
			tmp_pomcpow_iters = 100
			tmp_seed = 1
			k = (shape, (grid_dim,grid_dim,1), tmp_pomcpow_iters)
			config = results_regret[k][:config][tmp_seed]
			trial = initialize(convert(MEConfiguration, config))
			b0_set[k] = trial.b0
		end
	end
	b0_set
end

# ╔═╡ e6984433-c01e-4d55-b316-3b569c19dbb7
results_regret[(:blob, (10,10,1), 100)]

# ╔═╡ 5f358380-422f-4f0a-b59e-62337d08ebb9
# results_regret[(:blob, (10,10,1), 100)][:mass_map][4] |> MineralExploration.plot_ore_map

# ╔═╡ 7fc2d9ad-338b-438b-813e-1e9fc8f4ed72
results_regret[(:blob, (10,10,1), 100)][:rel_errors][1]

# ╔═╡ 295beb84-b6dd-4df0-aeed-623ac0ab6dfb
begin
	mean_distances = Float64[]
	mean_biases = Float64[]
	sterr_distances = Float64[]
	sterr_biases = Float64[]
	centers = []
	rkeys = []
	num_seeds = 500 # 500
	for shape in [:blob]
		for grid_dim in [10,30,50]
			ratio = (grid_dim,grid_dim,1) ./ (50,50,1)
			for pomcpow_iters in [100, 1_000, 10_000]
				b0_k = (shape, (grid_dim, grid_dim, 1), 100)
				k = (shape, (grid_dim, grid_dim, 1), pomcpow_iters)
				b0 = b0_set[b0_k]
				b_center = belief_center(b0)
				push!(centers, b_center)
				local_distances = []
				local_biases = []
				for seed in 1:num_seeds
					true_mass_map = results_regret[k][:mass_map][seed]
					distance = mean_distance_to_center(b_center, true_mass_map, ratio)
					bias = last(results_regret[k][:rel_errors][seed])
					# push!(local_distances, distance)
					# push!(local_biases, bias)
					push!(mean_distances, distance)
					push!(mean_biases, bias)
					push!(rkeys, (k,seed))
				end
				# push!(mean_distances, mean(local_distances))
				# push!(mean_biases, mean(local_biases))
				# push!(sterr_distances, std(local_distances)/sqrt(num_seeds))
				# push!(sterr_biases, std(local_biases)/sqrt(num_seeds))
			end
		end
	end
end

# ╔═╡ 4c426d81-dece-4c20-893c-5b5ef74fbe1c
d_max_idx = argmax(mean_distances)

# ╔═╡ 7b9646f1-557b-4dde-8e91-56c63c5317b8
mean_distances[d_max_idx]

# ╔═╡ d5f5b2c0-a6c8-4f71-8c4d-1b19995d7e36
mean_biases[d_max_idx]

# ╔═╡ cec4eaca-4d27-4df6-bfd2-37118159c44d
rk = rkeys[d_max_idx]

# ╔═╡ e5e9b25d-2ba2-436d-a8fb-a4003def3429
results_regret[rk[1]][:mass_map][rk[2]] |> MineralExploration.plot_ore_map

# ╔═╡ 6a9ee444-c73e-46d9-8a7a-1940335d1ed4
results_regret[rk[1]][:n_drills][rk[2]]

# ╔═╡ af2c880c-7cef-435b-bfa9-c6db3b340c82
histogram2d(mean_distances, mean_biases, show_empty=true)

# ╔═╡ 7db9ebf7-26a4-4606-ba2d-31773dced294
Xs = map(x->x[1], centers)

# ╔═╡ cf666432-7bad-4a9a-8804-615b09e893ec
function linear_regression(X, y)
    𝐗 = [ones(size(y)) X]
    θ = 𝐗\y
    return x -> [1;x]'θ
end

# ╔═╡ 303cbcfd-ee37-48e8-8910-0db9472cad78
lrf = linear_regression(mean_distances, mean_biases)

# ╔═╡ e61c0615-0551-421a-a5a8-19eec34744f1
begin
	scatter(mean_distances, mean_biases,
		xlabel="distance", ylabel="bias", label=false)
	plot!(xlims()[1]:xlims()[2], lrf, lw=2, c=:red, label=false)
end

# ╔═╡ d08c337c-5927-4bcc-a3da-fec5be06f2c7
scatter(mean_distances, mean_biases,
	xerr=sterr_distances, yerr=sterr_biases,
	xlabel="distance", ylabel="bias", label=false)

# ╔═╡ 8bc60e56-4d27-47e0-ab46-3d2b875404aa
md"""
# Random policy
"""

# ╔═╡ 80b89c42-016d-491d-a087-572b78671c0e
results_random = BSON.load("..\\scripts\\MEParallel.jl\\results\\results_random_policy.bson")[:results]

# ╔═╡ ce9ec1ff-6cb2-43cc-874d-ceb62bf083ae
begin
	random50x50blob = results_random[(:blob, (50,50,1), -1)]
	random50x50ellipse = results_random[(:ellipse, (50,50,1), -1)]
	random50x50circle = results_random[(:circle, (50,50,1), -1)]
end;

# ╔═╡ 98fdfbb9-0c0e-49aa-bac9-b98b4ae040f3
plot_rel_error_aggregate(random50x50blob, random50x50ellipse, random50x50circle)

# ╔═╡ 90a8ccd0-5d19-4105-973b-f23774250062
md"""
# Random policy (100 bores)
"""

# ╔═╡ c35ad770-9e03-4781-91a5-4676fcd623a1
results_random_100 = BSON.load("..\\scripts\\MEParallel.jl\\results\\results_random_policy_100.bson")[:results]

# ╔═╡ 7731ca9a-3bd5-4f7f-ac08-e43d0e62e52f
begin
	random50x50blob_100 = results_random_100[(:blob, (50,50,1), -1)]
	random50x50ellipse_100 = results_random_100[(:ellipse, (50,50,1), -1)]
	random50x50circle_100 = results_random_100[(:circle, (50,50,1), -1)]
end;

# ╔═╡ 9a751b96-7ca3-4aff-9df6-095c0538f6b5
plot_rel_error_aggregate(random50x50blob_100, random50x50ellipse_100, random50x50circle_100)

# ╔═╡ 1076dfdd-5bdf-48a1-8d30-04479ed5507e
random50x50blob_100

# ╔═╡ 766e2f39-cca7-4a74-8d5c-a1b327151e97
md"""
# Particle filter analysis
"""

# ╔═╡ 74623ee4-2c9f-46a2-b3d6-5a2a027251c7
up = MineralExploration.MEBeliefUpdater(trial.problem, 1000, 2.0) # trial.updater

# ╔═╡ 56b4340c-41c8-4544-9df7-3b9225ad173c
b = trial.b0

# ╔═╡ fb255509-357b-4d2c-8781-f106ade60f08
a = MineralExploration.MEAction(coords=CartesianIndex(10,18))

# ╔═╡ ca3da8f0-f7e2-43d6-afe1-d86be452a53f
	MineralExploration.plot_mass_map(trial.s0.ore_map, 0.7)[1]

# ╔═╡ f4992c17-6335-4303-a91a-f74c78ae8368
begin
	ore_quality = MineralExploration.high_fidelity_obs(trial.problem, trial.s0.ore_map, a)
	o = MineralExploration.MEObservation(ore_quality, false, false)
end

# ╔═╡ 5b48ebdc-0028-4d33-b913-d4e154e55ffd
MineralExploration.plot(b)

# ╔═╡ ec7b4eae-9ae1-42ea-9bc3-567ba8ab26fd
b2 = MineralExploration.update(up, b, a, o);

# ╔═╡ 7e438683-2a58-439e-bddd-b81e3385352e
MineralExploration.plot(b2)

# ╔═╡ 0f13ff64-077b-402f-9791-0a4b413529de
md"""
---
"""

# ╔═╡ 14a7334e-b2a8-4a45-b260-3536ab94a375
a2 = MineralExploration.MEAction(coords=CartesianIndex(10,15))

# ╔═╡ be39884d-ae63-470d-8e8e-e82ed1f3603d
begin
	MineralExploration.plot_ore_map(trial.s0.ore_map)
	a_hf = a.coords.I ./ ratio[1:2]
	scatter!([a_hf[2]], [a_hf[1]], label=false)
	a2_hf = a2.coords.I ./ ratio[1:2]
	scatter!([a2_hf[2]], [a2_hf[1]], label=false)
end

# ╔═╡ 2f3972e6-9478-452e-a429-7543090e6d5e
o2 = MineralExploration.MEObservation(
		MineralExploration.high_fidelity_obs(trial.problem, trial.s0.ore_map, a2),
	    false, false)

# ╔═╡ 237449d2-51bc-46e8-9ee5-8078465ed766
b3 = MineralExploration.update(up, b2, a2, o2)

# ╔═╡ cf4e6a86-b69d-451b-9fb2-c19fb954b73d
MineralExploration.plot(b3)

# ╔═╡ 85ab80d9-d9f4-4fe2-b952-1670fff3dd55
begin
	_, r_massive = MineralExploration.plot_mass_map(trial.s0.ore_map, trial.problem.massive_threshold)
	MineralExploration.plot_volume(trial.problem, b3, r_massive)[1]
end

# ╔═╡ 46435d91-cb3b-4b87-a00d-1712afe8e305
b3.particles[1].mainbody_params

# ╔═╡ e937fc37-e072-4156-afe1-3eb664816dfc
md"""
# Perturb
"""

# ╔═╡ 3685338c-790d-4deb-99a1-0bedc9cea12d
function create_gif(input_frames, output; fps=3)
	frames = Frames(MIME("image/png"), fps=fps)
	for frame in input_frames
		push!(frames, frame)
	end
	write(output, frames)
	LocalResource("./$output")
end

# ╔═╡ fe444c99-c62b-4df3-9666-2f0227f0ffbb
trial.s0.mainbody_map |> MineralExploration.plot_ore_map

# ╔═╡ e9723f21-efd1-461f-8bbf-398572f795af
trial.s0.mainbody_params

# ╔═╡ 345b123b-835f-42ba-bace-32274f146877
begin
	frames = []
	global p_mb
	global gif_prefix
	use_blob = true
	for i in 1:20
		global p_mb, gif_prefix

		if use_blob
			mb = BlobNode(grid_dims=(50,50,1),
				# center=trial.s0.mainbody_params[1],
				# N=trial.s0.mainbody_params[2],
				# factor=trial.s0.mainbody_params[3],
				# points=trial.s0.mainbody_params[4],
				# angle=trial.s0.mainbody_params[5]
				)
			p_mb, p_mb_params = 
				MineralExploration.perturb_sample(mb, trial.s0.mainbody_params, 2; 
					recompute_points=false, copy_points=false, clamped_to_prior=true)
			title = "blob perturbation"
			gif_prefix = "blob"
		else
			mb = EllipseNode(grid_dims=(50,50,1))
			_, mb_params = rand(mb)
			p_mb, p_mb_params = 
				MineralExploration.perturb_sample(mb, mb_params, 2)
			title = "ellipse perturbation"
			gif_prefix = "ellipse"
		end
		frame = MineralExploration.plot_ore_map(p_mb)
		title!(title)
		push!(frames, frame)
	end
end

# ╔═╡ ba1d0d6a-50e7-4900-be31-9d1a151241e3
p_mb |> MineralExploration.plot_ore_map

# ╔═╡ f5700ac2-4646-446e-8e1e-106898e51cbc
create_gif(frames, "$(gif_prefix)_perturbed_clamped.gif")

# ╔═╡ Cell order:
# ╠═a333c86e-d90e-49e5-8623-8e844ad2a9b8
# ╟─fb8868be-997f-11ec-0fb4-37f4c70c68c7
# ╠═62013d0f-8cbd-473a-8ed8-9662a567328b
# ╠═7d875d7d-e2e8-413f-aece-452b539fca16
# ╠═ba7a7f30-fc84-4a5e-b853-9946c69a1920
# ╠═d1c80a9e-e5a9-458a-b76f-8dfaa5ab8b53
# ╠═fe64f46b-98f6-454a-9ab8-2f73a9722e11
# ╠═3ac672e1-33de-4405-ace2-b982ab8f94df
# ╟─4185c93a-42b1-43b8-b44d-7c301b270392
# ╠═77ca5cd5-1501-48d9-a9ba-f3b0c482c8f2
# ╟─06a55c95-92cc-453f-8e10-1d53f153fb23
# ╠═4a144a76-625c-4b17-b998-89b7fef5014d
# ╟─80bbd217-d81e-4b80-b455-a84309c4fe8d
# ╠═fc372d8f-2f97-4f5c-8062-c35f7864912f
# ╟─fa227921-0730-4f33-ae61-8ea635be9e44
# ╠═ed628d8b-cf6e-4e4a-865b-d844f676d512
# ╠═916b42b0-b171-44fc-b27c-960208c2f843
# ╠═4b00b8cd-b58c-40a1-b86b-adf77918ec3c
# ╠═ff7463c2-1f8f-44ae-85e0-fff247bd5560
# ╟─7a23a34c-a0e4-47dd-970f-f42474bde9c5
# ╠═321013a9-6371-4628-adaf-ca89017041e0
# ╠═9f32913d-7913-47f5-90f6-dbfa6955e46d
# ╠═1ce9a2c4-917e-47e7-bb0b-4000aed7ae3d
# ╠═dd1fab4d-1c10-4a55-942b-45cdc3a3b96f
# ╠═34d07754-d716-4c4d-a346-56b73947b965
# ╠═56f2a8ec-192a-4d5b-ae82-376062bbac8e
# ╠═d4a91997-9e5d-4b3e-809d-f0532f06e07d
# ╠═ebb31bbf-cee0-4a09-a044-b864e0856e84
# ╠═816a3b19-8e43-4b84-86fc-893c4ed92619
# ╠═3a33df07-bcff-4f22-abb8-565613cd8fa9
# ╠═74555820-8c6e-4714-86f2-d6c33434bb80
# ╠═729baffb-154b-4e30-8918-8547c1071177
# ╠═9c1c3ebf-a63a-41c1-9f7d-6d853e282a74
# ╠═776b4f28-81ff-44ef-b6c2-8b418075cc1e
# ╠═88c4fc71-fb7d-4b44-9cfd-77c100437c3a
# ╠═1d4e6d83-e96b-4138-b6c1-844b20fe3729
# ╠═e74da8ce-ff04-4778-9ec6-25ce30ed7cd2
# ╠═f253b3b3-bcac-4330-af1a-6cd098d4375f
# ╠═6b53029d-ffe8-4fb8-ac97-9b1b5677821e
# ╠═7bc1f9c8-b771-4538-8970-05c4d4f5d12f
# ╠═4f6d7b0a-e424-4348-a566-5300a0ce1420
# ╠═44a5c6ca-59de-4577-82c5-8fb87cdb6a47
# ╠═296c6deb-1ce4-4c9e-af85-2977ffa29990
# ╠═8daff9ab-c4fc-4d34-95fc-a47ab70aec50
# ╠═b15ad0b4-37c9-4533-9d1c-5bc797079043
# ╟─71782ce8-1fdf-4140-a53a-9ad70ce5779f
# ╠═506dc301-81e8-4c80-bf7e-3843f9b911d2
# ╠═6acd1e62-ca92-4f62-96a6-38bc71ec6871
# ╠═6b756e10-0f71-4ccb-bc23-ff9a74bba135
# ╠═2f03b56f-4ea7-4dd3-8974-ecd570ff6c0b
# ╠═2a7d356e-20f8-403b-9e9b-726e54919643
# ╠═b850b7fc-641e-498d-87ff-43f6e938abb6
# ╠═d28a33db-b17c-415d-923b-8cdf87f5c860
# ╠═de18f869-8f5e-4c75-b631-da8ae0f61794
# ╠═eb9cd498-27b6-4747-9d84-5b02b73c9eae
# ╠═c781d26c-bbcc-4ea9-a2d4-128d3dfd9bd4
# ╠═235d6716-2908-497f-9001-4b1ed65b5a81
# ╠═9d5c8911-3934-465f-9337-b12ef66ca7b2
# ╠═6387dd41-dfc4-4109-8602-077dfb13757e
# ╟─2a4125ad-2e2a-4e12-ac38-04ed56c45c2c
# ╠═081c6ba0-4155-4d96-a8d8-b387f2cc5244
# ╠═a45ba4d6-df96-47c5-a129-a6eab357bed3
# ╟─ae7606ef-fc61-403d-8136-9e1ce23a6a84
# ╠═deddcc30-62e4-4736-830b-550d880970d1
# ╠═28840030-6459-4c6f-a38e-470eb5d9c1c6
# ╠═fbf32114-5237-4997-96c4-cc5475862bdd
# ╠═56b40041-50fa-43f3-b260-f3b189bf21a7
# ╠═211cbb7d-d2d8-4e5b-9fbe-6a90bc543985
# ╠═4c1554ea-7451-4edd-b2d6-d4c56d17d7f0
# ╟─675f654e-ad55-432a-a421-ca9a28737353
# ╠═e11f7f67-ac9b-4df5-a678-8172b499152f
# ╟─c96dd779-621a-4c16-b4eb-375b4d67c299
# ╠═4bec3cad-4e8f-4b59-b878-b584cf8daacd
# ╠═763e4c58-72a6-4286-b9ea-c825309105dd
# ╠═b8d17852-8668-42ca-98fa-e7f1424ffeb2
# ╠═55433d68-5dd6-4d1e-a0c7-abde02d7cf65
# ╠═61106a23-8901-450f-a58a-05c58163f3a8
# ╠═b8450d7d-ab0f-4429-9a82-169d4eda5fd7
# ╠═db0123bc-103e-4ab5-a60f-02bef918dfd9
# ╠═19ecdfb4-4cfa-445b-af60-92d376cd7fd1
# ╠═6cd5b0cd-685c-4c89-af98-59664a64c95b
# ╠═b8a6cf63-04ee-480f-82a4-f2f8a9cb20b1
# ╠═8117557b-035d-49bc-801a-11cff86208e3
# ╠═3ba758ab-d8da-4f67-9b46-b941fe35f4fe
# ╠═2eb69b52-bc1b-469d-8829-de6108d95a01
# ╠═1c6b50e0-7024-45a2-aa6b-02031eb880a4
# ╠═da2e51af-1b80-4927-ba7c-b361891662f7
# ╠═2576c794-419b-4de0-8b71-6f5c50a2b9cb
# ╠═80f8f9ca-1edb-4305-9aac-178bd494c078
# ╠═1bf3c03a-34e7-4a8e-934a-691afa0586d6
# ╠═ed2e6025-aa46-4565-8703-e5ad4a8656ce
# ╠═d32a110e-ca3e-4ba2-8d92-baf93883f67a
# ╠═3ec60a22-40b4-4750-88cb-509621aa1c52
# ╠═303198db-a517-483f-abc1-ff27f921636a
# ╠═b0258367-59ce-4b98-b506-f2de7a64f333
# ╠═fe849f40-ace8-41f7-9f10-099209342d67
# ╠═7d41f3be-20d8-4726-b0f9-fd732e4e4abe
# ╠═f474bc29-4161-4e66-8ab5-a41f3cdaf670
# ╠═f78a4a1a-8075-4393-8736-5f5dfd61c237
# ╠═729e1916-92c6-447c-bfd0-8094820feff8
# ╠═41ae1023-32c1-4dac-b134-2b01ddf2863d
# ╠═f77a8482-e4d6-47ab-b186-98a2e3a7ed03
# ╠═8a1d7ee7-95db-49bd-a9e2-51f17ea07a6a
# ╠═c2fa1f0f-4ba3-4269-9ff4-be84efc06424
# ╠═ff9561f5-ead3-4849-a97e-805ba844d60a
# ╠═e6984433-c01e-4d55-b316-3b569c19dbb7
# ╠═5f358380-422f-4f0a-b59e-62337d08ebb9
# ╠═7fc2d9ad-338b-438b-813e-1e9fc8f4ed72
# ╠═295beb84-b6dd-4df0-aeed-623ac0ab6dfb
# ╠═4c426d81-dece-4c20-893c-5b5ef74fbe1c
# ╠═7b9646f1-557b-4dde-8e91-56c63c5317b8
# ╠═d5f5b2c0-a6c8-4f71-8c4d-1b19995d7e36
# ╠═cec4eaca-4d27-4df6-bfd2-37118159c44d
# ╠═e5e9b25d-2ba2-436d-a8fb-a4003def3429
# ╠═6a9ee444-c73e-46d9-8a7a-1940335d1ed4
# ╠═e61c0615-0551-421a-a5a8-19eec34744f1
# ╠═af2c880c-7cef-435b-bfa9-c6db3b340c82
# ╠═7db9ebf7-26a4-4606-ba2d-31773dced294
# ╠═cf666432-7bad-4a9a-8804-615b09e893ec
# ╠═303cbcfd-ee37-48e8-8910-0db9472cad78
# ╠═d08c337c-5927-4bcc-a3da-fec5be06f2c7
# ╟─8bc60e56-4d27-47e0-ab46-3d2b875404aa
# ╠═80b89c42-016d-491d-a087-572b78671c0e
# ╠═ce9ec1ff-6cb2-43cc-874d-ceb62bf083ae
# ╠═98fdfbb9-0c0e-49aa-bac9-b98b4ae040f3
# ╟─90a8ccd0-5d19-4105-973b-f23774250062
# ╠═c35ad770-9e03-4781-91a5-4676fcd623a1
# ╠═7731ca9a-3bd5-4f7f-ac08-e43d0e62e52f
# ╠═9a751b96-7ca3-4aff-9df6-095c0538f6b5
# ╠═1076dfdd-5bdf-48a1-8d30-04479ed5507e
# ╟─766e2f39-cca7-4a74-8d5c-a1b327151e97
# ╠═74623ee4-2c9f-46a2-b3d6-5a2a027251c7
# ╠═56b4340c-41c8-4544-9df7-3b9225ad173c
# ╠═fb255509-357b-4d2c-8781-f106ade60f08
# ╠═be39884d-ae63-470d-8e8e-e82ed1f3603d
# ╠═ca3da8f0-f7e2-43d6-afe1-d86be452a53f
# ╠═f4992c17-6335-4303-a91a-f74c78ae8368
# ╠═5b48ebdc-0028-4d33-b913-d4e154e55ffd
# ╠═ec7b4eae-9ae1-42ea-9bc3-567ba8ab26fd
# ╠═7e438683-2a58-439e-bddd-b81e3385352e
# ╟─0f13ff64-077b-402f-9791-0a4b413529de
# ╠═14a7334e-b2a8-4a45-b260-3536ab94a375
# ╠═2f3972e6-9478-452e-a429-7543090e6d5e
# ╠═237449d2-51bc-46e8-9ee5-8078465ed766
# ╠═cf4e6a86-b69d-451b-9fb2-c19fb954b73d
# ╠═85ab80d9-d9f4-4fe2-b952-1670fff3dd55
# ╠═46435d91-cb3b-4b87-a00d-1712afe8e305
# ╟─e937fc37-e072-4156-afe1-3eb664816dfc
# ╠═3685338c-790d-4deb-99a1-0bedc9cea12d
# ╠═fe444c99-c62b-4df3-9666-2f0227f0ffbb
# ╠═e9723f21-efd1-461f-8bbf-398572f795af
# ╠═345b123b-835f-42ba-bace-32274f146877
# ╠═ba1d0d6a-50e7-4900-be31-9d1a151241e3
# ╠═f5700ac2-4646-446e-8e1e-106898e51cbc
