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

# ╔═╡ 4d723203-f873-4ce9-a784-7ebaeca42f73
begin
	using Revise
	using Random
	using Pkg
	Pkg.develop(path="../../MineralExploration")
	using MineralExploration
end

# ╔═╡ 18aedc2e-73aa-4bc5-91f9-dacc74033587
using Plots; default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

# ╔═╡ f62f97d3-8493-45b5-9fee-c908a57a6178
using POMDPs

# ╔═╡ 53c56583-4f56-4046-8b60-3d82a9d626bb
using Images

# ╔═╡ c2794b75-75b2-42ce-8db0-4706d175be7e
using Parameters, Distributions

# ╔═╡ 18ba3d34-e835-4c97-a986-21253a5dc57b
using PlutoUI

# ╔═╡ 3d607e3c-1af4-4ea1-b443-acb4dc0df752
md"""
# Mainbody fidelity
"""

# ╔═╡ 22aca29a-c41b-409b-a85b-ae5cb06ece16
begin
	N_INITIAL = 0
	MAX_BORES = 25
	MIN_BORES = 10
	GRID_SPACING = 1
	
	hf_grid_dims = (50, 50, 1)
	grid_dims = (80, 80, 1)
	ratio = grid_dims ./ hf_grid_dims
	hf_mainbody_loc_bounds = [12.0, 50.0]
	mainbody_loc_bounds = hf_mainbody_loc_bounds .* ratio[1:2]
	
	MAX_MOVEMENT = 10
	mainbody = MultiVarNode(grid_dims=grid_dims,
						    mainbody_loc_bounds=mainbody_loc_bounds)
end

# ╔═╡ 9d3a7c27-465b-48e7-90a9-c989d59f7eca
begin
	m = MineralExplorationPOMDP(
		max_bores=MAX_BORES,
		delta=GRID_SPACING+1,
		grid_spacing=GRID_SPACING,
		mainbody_gen=mainbody,
		max_movement=MAX_MOVEMENT,
		min_bores=MIN_BORES,
		grid_dim=grid_dims)

	initialize_data!(m, N_INITIAL)
end

# ╔═╡ 7564032c-1154-4052-aea9-14b471db92eb
function plot_map(A; clims=true)
	xl = (1,size(A,1))
	yl = (1,size(A,2))	
	cl = clims ? (0,1) : (minimum(A), maximum(A))
	heatmap(A, fill=true, clims=cl, aspect_ratio=1, xlims=xl, ylims=yl)
end

# ╔═╡ 0c077c46-dbe4-4de2-8064-1edf8cfef8d1
function bitmap(A, threshold=0.3)
	return BitMatrix(A .> threshold)
end

# ╔═╡ e10f53bf-f1c1-4110-9c51-c384066052c2
## Multiple Variable Node (Elliptic)
@with_kw struct MultiVarNodeElliptic <: MainbodyGen
    grid_dims::Tuple{Int64, Int64, Int64} = (50, 50, 1)
    mainbody_loc_bounds::Vector{Float64} = [15.0, 35.0]
    mainbody_var_min::Float64 = 20.0
    mainbody_var_max::Float64 = 40.0
    n_nodes::Int64 = 2
	stretch::Vector{Vector{Float64}} = [[1, 4], [4, 1]]
end

# ╔═╡ 0199dde0-dc37-42be-bfc6-294c73a696f7
function Base.rand(rng::Random.AbstractRNG, mb::MultiVarNodeElliptic)
    x_dim = mb.grid_dims[1]
    y_dim = mb.grid_dims[2]
    lode_map = zeros(Float64, x_dim, y_dim)
    mainbody_params = []
    mainbody_var = rand(rng)*(mb.mainbody_var_max - mb.mainbody_var_min) + mb.mainbody_var_min
    for idx in 1:mb.n_nodes
        d_lode_map = zeros(Float64, x_dim, y_dim)
        rand_loc = rand(2)
        cov = Distributions.PDiagMat(mb.stretch[idx] .* [mainbody_var,mainbody_var])
		mainbody_loc = rand_loc .* (mb.mainbody_loc_bounds[2] - mb.mainbody_loc_bounds[1]) .+ mb.mainbody_loc_bounds[1]
		mvnorm = MvNormal(mainbody_loc, cov)
        for i = 1:x_dim
            for j = 1:y_dim
                d_lode_map[i, j] = pdf(mvnorm, [float(i), float(j)])
            end
        end
        d_mb_params = vcat(mainbody_loc, mainbody_var)
        lode_map += d_lode_map
        push!(mainbody_params, d_mb_params)
    end
    # lode_map ./= mb.n_nodes
    return (lode_map, mainbody_params)
end

# ╔═╡ 4c6ea470-f32d-4a99-850d-2bc4ed50b0d3
Base.rand(mb::MultiVarNodeElliptic) = rand(Random.GLOBAL_RNG, mb)

# ╔═╡ e77dde66-50f8-4b1d-90bf-108f51b4e890
begin
	Random.seed!(0)
	ds0 = POMDPs.initialstate_distribution(m)
	s0 = rand(ds0)
end;

# ╔═╡ ffb78c8c-63bb-4656-bc60-c8e48f82630e
ore_map = s0.ore_map[:,:,1];

# ╔═╡ d42049c0-0ee2-4247-9f93-39404a80f0f6
plot_map(ore_map); title!("ore map")

# ╔═╡ 8a3fa7e4-acc4-4ccd-8e55-8f115a5d0300
ore_map_resize = imresize(ore_map, (10, 10))

# ╔═╡ 3d5944f4-f9be-4dbf-aa6b-ef920e410099
plot_map(ore_map_resize)

# ╔═╡ 82a757e8-a9a7-4ca4-bf36-60bc09bb6f7d
mainbody_map = s0.mainbody_map[:,:,1]

# ╔═╡ 97347c80-fee0-474e-9fee-b3d80d7e17f4
plot_map(mainbody_map); title!("mainbody map")

# ╔═╡ 82abcf3e-d02f-4a41-8ac4-2afe8d77636b
plot_map(mainbody_map)

# ╔═╡ 350e511e-9e9f-49a0-be10-f43995afd829
A = mainbody_map

# ╔═╡ 691b6c45-8c83-4424-846c-d416acc04915
bitmap(mainbody_map) |> plot_map

# ╔═╡ 95ca522b-0f92-4b0a-8c04-f60f3bdcdf17
bitmap(imresize(bitmap(mainbody_map), (10,10))) |> plot_map

# ╔═╡ 5697d72e-eafb-4b96-87c0-b35eb368a1d8
mb = MultiVarNodeElliptic(grid_dims=grid_dims,
	mainbody_loc_bounds=mainbody_loc_bounds)

# ╔═╡ 4f0044bb-124b-48a5-83d0-e05d19e05995
rand(mb)[1] |> A->heatmap(A, aspect_ratio=1, xlims=(1, grid_dims[1]))

# ╔═╡ 4cf46d25-0185-4ae3-839c-a5ca4adee5a4
Distributions.PDiagMat([28, 28])

# ╔═╡ 3ac38dce-2820-4a7d-b572-0e5b0dba6620
md"""
## Mainbody (MvNormal)
"""

# ╔═╡ 49ebc864-0d72-4f32-9b86-2cb0ca0b8ae6
mainbody

# ╔═╡ 2806990f-6c3b-489a-b0bb-45e2d719f904
lode_map, mainbody_var = rand(mainbody)

# ╔═╡ 89c87e2f-7522-4ecf-abe3-93105eed28a4
@bind is_perturbed CheckBox()

# ╔═╡ 1bddd84d-98dd-4333-bb25-fa3d704262ae
p_lode_map, p_mainbody_params = MineralExploration.perturb_sample(mainbody, mainbody_var, 1.0)

# ╔═╡ 2b96b92c-f0ef-45f6-96d8-1547893d0e95
if is_perturbed
	plot_map(p_lode_map, clims=false)
	title!("perturbed")
else
	plot_map(lode_map, clims=false)
	title!("original")
end

# ╔═╡ 0e01c62f-25c1-4794-bb73-76ee41d1ee0e
mainbody_var .- p_mainbody_params

# ╔═╡ 6d8f32a3-ab86-48d0-8649-bfdf194cd577
md"""
# Single
"""

# ╔═╡ 4019d210-3038-4752-a7b1-d55116792a1c
single_mainbody = SingleVarNode(grid_dims=grid_dims)

# ╔═╡ e559e2ce-b0b4-41c2-88a4-6238fc8e894f
single_mainbody_map, single_mainbody_var = rand(single_mainbody)

# ╔═╡ 8587401d-f05e-4d52-86fc-a25e64770fc5
plot_map(single_mainbody_map, clims=false)

# ╔═╡ 8854be14-5db6-44e2-9427-42bd6e3154db
single_mainbody_var

# ╔═╡ Cell order:
# ╟─3d607e3c-1af4-4ea1-b443-acb4dc0df752
# ╠═18aedc2e-73aa-4bc5-91f9-dacc74033587
# ╠═4d723203-f873-4ce9-a784-7ebaeca42f73
# ╠═f62f97d3-8493-45b5-9fee-c908a57a6178
# ╠═22aca29a-c41b-409b-a85b-ae5cb06ece16
# ╠═9d3a7c27-465b-48e7-90a9-c989d59f7eca
# ╠═e77dde66-50f8-4b1d-90bf-108f51b4e890
# ╠═7564032c-1154-4052-aea9-14b471db92eb
# ╠═ffb78c8c-63bb-4656-bc60-c8e48f82630e
# ╠═d42049c0-0ee2-4247-9f93-39404a80f0f6
# ╠═82a757e8-a9a7-4ca4-bf36-60bc09bb6f7d
# ╠═97347c80-fee0-474e-9fee-b3d80d7e17f4
# ╠═53c56583-4f56-4046-8b60-3d82a9d626bb
# ╠═8a3fa7e4-acc4-4ccd-8e55-8f115a5d0300
# ╠═3d5944f4-f9be-4dbf-aa6b-ef920e410099
# ╠═82abcf3e-d02f-4a41-8ac4-2afe8d77636b
# ╠═350e511e-9e9f-49a0-be10-f43995afd829
# ╠═691b6c45-8c83-4424-846c-d416acc04915
# ╠═95ca522b-0f92-4b0a-8c04-f60f3bdcdf17
# ╠═0c077c46-dbe4-4de2-8064-1edf8cfef8d1
# ╠═c2794b75-75b2-42ce-8db0-4706d175be7e
# ╠═e10f53bf-f1c1-4110-9c51-c384066052c2
# ╠═0199dde0-dc37-42be-bfc6-294c73a696f7
# ╠═4c6ea470-f32d-4a99-850d-2bc4ed50b0d3
# ╠═5697d72e-eafb-4b96-87c0-b35eb368a1d8
# ╠═4f0044bb-124b-48a5-83d0-e05d19e05995
# ╠═4cf46d25-0185-4ae3-839c-a5ca4adee5a4
# ╟─3ac38dce-2820-4a7d-b572-0e5b0dba6620
# ╠═18ba3d34-e835-4c97-a986-21253a5dc57b
# ╠═49ebc864-0d72-4f32-9b86-2cb0ca0b8ae6
# ╠═2806990f-6c3b-489a-b0bb-45e2d719f904
# ╠═89c87e2f-7522-4ecf-abe3-93105eed28a4
# ╠═2b96b92c-f0ef-45f6-96d8-1547893d0e95
# ╠═1bddd84d-98dd-4333-bb25-fa3d704262ae
# ╠═0e01c62f-25c1-4794-bb73-76ee41d1ee0e
# ╟─6d8f32a3-ab86-48d0-8649-bfdf194cd577
# ╠═4019d210-3038-4752-a7b1-d55116792a1c
# ╠═e559e2ce-b0b4-41c2-88a4-6238fc8e894f
# ╠═8587401d-f05e-4d52-86fc-a25e64770fc5
# ╠═8854be14-5db6-44e2-9427-42bd6e3154db
