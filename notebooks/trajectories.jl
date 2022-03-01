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

# ╔═╡ deddcc30-62e4-4736-830b-550d880970d1
using Distributions

# ╔═╡ fb8868be-997f-11ec-0fb4-37f4c70c68c7
md"""
# Trajectories
"""

# ╔═╡ 62013d0f-8cbd-473a-8ed8-9662a567328b
results_regret = BSON.load("..\\scripts\\MEParallel.jl\\results\\results_regret.bson")[:results]

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
τ_min_accuracy = fargmin(results_regret, accuracy)

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

# ╔═╡ 71782ce8-1fdf-4140-a53a-9ad70ce5779f
md"""
# Worst-case: Accuaracy
"""

# ╔═╡ 506dc301-81e8-4c80-bf7e-3843f9b911d2
τ_min_accuracy

# ╔═╡ 6acd1e62-ca92-4f62-96a6-38bc71ec6871
wc_acc_res = results_regret[τ_min_accuracy]

# ╔═╡ 2f03b56f-4ea7-4dd3-8974-ecd570ff6c0b
wc_acc_idx = argmax(regret(results_regret[τ_min_accuracy]))

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

# ╔═╡ de18f869-8f5e-4c75-b631-da8ae0f61794
config = convert(MEConfiguration, wc_acc_config)

# ╔═╡ eb9cd498-27b6-4747-9d84-5b02b73c9eae
trial = initialize(config)

# ╔═╡ c781d26c-bbcc-4ea9-a2d4-128d3dfd9bd4
results = evaluate(trial; save_dir="worst_case_accuracy")

# ╔═╡ 6387dd41-dfc4-4109-8602-077dfb13757e
results.ore_map |> MineralExploration.plot_ore_map

# ╔═╡ 2ba1013f-7262-42ac-b07a-a9d128f1ce17
wc_acc_res

# ╔═╡ 235d6716-2908-497f-9001-4b1ed65b5a81
results_dict = MEParallel.convert(Dict, results)

# ╔═╡ 9d5c8911-3934-465f-9337-b12ef66ca7b2
regret(results_dict)

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
	histogram(regret(results[k]), bins=4, normalize=:pdf, c=:gray)
	𝒩 = fit(Normal, regret(results[k]))
	𝒩 = TruncatedNormal(𝒩.μ, 𝒩.σ, 0, Inf)
	plot!(x->pdf(𝒩,x), c=:red, lw=2)
	title!(string(k))
end

# ╔═╡ 4c1554ea-7451-4edd-b2d6-d4c56d17d7f0
regret_distribution(results_regret, τ_max_expected_regret)

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
# ╟─71782ce8-1fdf-4140-a53a-9ad70ce5779f
# ╠═506dc301-81e8-4c80-bf7e-3843f9b911d2
# ╠═6acd1e62-ca92-4f62-96a6-38bc71ec6871
# ╠═2f03b56f-4ea7-4dd3-8974-ecd570ff6c0b
# ╠═2a7d356e-20f8-403b-9e9b-726e54919643
# ╠═b850b7fc-641e-498d-87ff-43f6e938abb6
# ╠═d28a33db-b17c-415d-923b-8cdf87f5c860
# ╠═de18f869-8f5e-4c75-b631-da8ae0f61794
# ╠═eb9cd498-27b6-4747-9d84-5b02b73c9eae
# ╠═c781d26c-bbcc-4ea9-a2d4-128d3dfd9bd4
# ╠═6387dd41-dfc4-4109-8602-077dfb13757e
# ╠═2ba1013f-7262-42ac-b07a-a9d128f1ce17
# ╠═235d6716-2908-497f-9001-4b1ed65b5a81
# ╠═9d5c8911-3934-465f-9337-b12ef66ca7b2
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
