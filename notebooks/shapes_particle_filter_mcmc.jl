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

# ╔═╡ 7fb2c2da-ca38-489c-b8c1-1cc402b68222
begin
	using Revise
	using Random
end

# ╔═╡ 909a415e-f420-441c-bd2d-4ec86239113b
using Plots; default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

# ╔═╡ 7062f835-165c-414b-b25f-a060de10b3cc
using POMDPs

# ╔═╡ 58dab9b5-25b2-42fe-b918-0dea561d5fd5
using Images

# ╔═╡ a4aa5758-6caf-47eb-a18c-ff2281453e23
using Parameters

# ╔═╡ 4d70219e-9fd3-4d14-80d2-3c12fff9e837
using Distributions

# ╔═╡ 718b84cd-9851-4b4a-b9fb-b0d1ec2a628e
using Luxor

# ╔═╡ e78b494c-d32c-4123-a19c-9b1ec2fca78b
using ImageFiltering

# ╔═╡ e487916b-3581-4814-b2cd-36939fb471c7
using PlutoUI

# ╔═╡ 522efd58-3756-4ad7-931e-ef1f0bd8c483
using ProgressLogging

# ╔═╡ 2039b5c2-b834-495a-a3dd-ea20394f9bef
using LinearAlgebra

# ╔═╡ 668a9810-7a3f-11ec-0861-89fc669dae9e
md"""
# Mainbody fidelity
"""

# ╔═╡ bd1330a9-b65a-46a5-93ce-ecaf29632699
massive_threshold = 0.7

# ╔═╡ 809bb5c7-8635-451a-8bba-11c3b5da215b
md"""
## Shapes
Non-circular shapes includes a rotation angle $\alpha$
- Square: side length $\ell$
- Circle: radius $r$
- Ellipse: width $w$ and height $h$
"""

# ╔═╡ b216e135-e2a2-4646-9df4-1ccdda7d2dec
dims = (80,80)

# ╔═╡ 1ed349dc-a7d1-40b9-8c2b-e915d72fd150
# noise = imfilter(rand(dims...), Kernel.gaussian(1))

# ╔═╡ a91efc56-40e2-4a40-af75-b60d06f933c6
function limits(dims)
	lims_x = (0, dims[1]) .+ 1/2
	lims_y = (0, dims[2]) .+ 1/2
	return (lims_x, lims_y)
end

# ╔═╡ 28dd3755-9c37-4e30-b58f-4b04c3ded69e
lims_x, lims_y = limits(dims)

# ╔═╡ cee46dc3-32d4-4e79-96c4-5480f1402b80
# heatmap(noise, aspect_ratio=1, xlims=lims_x, ylims=lims_y, colorbar=false, c=:viridis)

# ╔═╡ 929486d5-1e98-4e03-a2d6-3279699d95df
md"""
# Luxor
"""

# ╔═╡ e26ba54c-efd3-44af-b8df-79569a31c160
abstract type AbstractShape end

# ╔═╡ ee6b7a25-2232-4130-81ab-5ddcf94e0b4f
abstract type Shape <: AbstractShape end

# ╔═╡ 355aebe8-a750-434e-8344-1d8b5ecc33fe
function plot_subsurface(noise, ore; clims=true, title="")
	dims = size(ore)
	lims_x, lims_y = limits(dims)
	subsurface = noise + ore
	heatmap(subsurface,
		aspect_ratio=1,
		xlims=lims_x,
		ylims=lims_y,
		c=:viridis,
		clims=clims ? (0,1) : (minimum(subsurface), maximum(subsurface)),
		title=title)
end

# ╔═╡ 8113c944-287f-442f-8c43-a00fb39b9ca5
function plot_subsurface(ore; kwargs...)
	return plot_subsurface(zeros(size(ore)...), ore; kwargs...)
end

# ╔═╡ a7d933d7-2e99-4303-bf8a-7b57981f5349
md"""
## Circles
"""

# ╔═╡ ce783848-56cc-43d6-907c-74dfd6d8c0cd
function center_distribution(dims, padding=0.2)
	xmargin = dims[1] * padding
	ymargin = dims[2] * padding
	return Product(Uniform.([1,1] .+ xmargin, [dims[1], dims[2]] .- ymargin))
end

# ╔═╡ f395c41f-399f-48f6-9b8c-3010d555887d
begin
	@with_kw struct Circle <: Shape
		dims::Tuple{Real, Real}
		center::Union{Distribution, Vector} = center_distribution(dims)
		radius::Union{Distribution, Real} = Uniform(dims[1]/80, dims[1]/8)
	end
	
	Circle(dims) = Circle(dims=dims)
end

# ╔═╡ 14e5bbf4-f51b-4d7a-9eba-cbe71b99fe61
Circle <: Shape <: AbstractShape

# ╔═╡ 3d8da47d-bee1-46e8-96ed-07e1c5ad897d
function draw(shape::Circle)
	return circle(Point(0,0), shape.radius, :fill)
end

# ╔═╡ 8f5e6e85-83a8-4410-8f8f-d84b6925b549
function Base.rand(shape::Circle)
	dims = shape.dims
	center = isa(shape.center, Distribution) ? rand(shape.center) : shape.center
	radius = isa(shape.radius, Distribution) ? rand(shape.radius) : shape.radius
	return Circle(dims=dims, center=center, radius=radius)
end

# ╔═╡ 3dd84717-7214-4d70-b82d-d1e0558c129a
md"""
## Ellipsis
"""

# ╔═╡ bc48337a-8699-45a5-a135-de75933e07a4
begin
	@with_kw struct Ellipse <: Shape
		dims::Tuple{Real, Real}
		center::Union{Distribution, Vector} = center_distribution(dims)
		width::Union{Distribution, Real} = Uniform(dims[1]/16, dims[1]/4)
		height::Union{Distribution, Real} = Uniform(dims[2]/16, dims[2]/4)
		angle::Union{Distribution, Real} = Uniform(0, 2π)
	end
	
	Ellipse(dims) = Ellipse(dims=dims)
end

# ╔═╡ 327409bc-bbf4-4236-9de5-19a61b69b9a1
function draw(shape::Ellipse)
	return ellipse(Point(0,0), shape.width, shape.height, :fill)
end

# ╔═╡ 06e993f6-4300-479b-8f4f-3c1666f29acc
function Base.rand(shape::Ellipse)
	dims = shape.dims
	center = isa(shape.center, Distribution) ? rand(shape.center) : shape.center
	width = isa(shape.width, Distribution) ? rand(shape.width) : shape.width
	height = isa(shape.height, Distribution) ? rand(shape.height) : shape.height
	angle = isa(shape.angle, Distribution) ? rand(shape.angle) : shape.angle
	return Ellipse(dims=dims, center=center, width=width, height=height, angle=angle)
end

# ╔═╡ 63799459-b343-4162-a421-ba3a77d50b65
seed_ellipse = 0

# ╔═╡ 8fc9b648-3810-4b32-8da8-7bf14225c6e8
md"""
## Bezier Blobs
"""

# ╔═╡ 2014a964-ab02-4fc4-962f-2f057a1da115
begin
	@with_kw struct Blob <: Shape
		dims::Tuple{Real, Real}
		center::Union{Distribution, Vector} = center_distribution(dims)
		N::Union{Distribution, Real} = Uniform(3, 10)
		factor::Union{Distribution, Real} = Uniform(3, 5)
		points::Union{Nothing, Vector{Point}} = nothing
		angle::Union{Distribution, Real} = Uniform(0, 2π)
	end

	Blob(dims) = Blob(dims=dims)
end

# ╔═╡ 9e7ca98c-34a2-44bb-ad72-e66e12dc006d
function draw(shape::Blob)
	dims = shape.dims
	factor = shape.factor
	pts = shape.points
	return drawbezierpath(makebezierpath(pts), :fill)
end

# ╔═╡ 336f0ea3-ea05-4207-bf2a-e21ac380f424
function Base.rand(shape::Blob)
	dims = shape.dims
	center = isa(shape.center, Distribution) ? rand(shape.center) : shape.center
	N = isa(shape.N, Distribution) ? rand(shape.N) : shape.N
	factor = isa(shape.factor, Distribution) ? rand(shape.factor) : shape.factor
	if isnothing(shape.points)
		pts = polysortbyangle(randompointarray(
			Point(0,0), Point(dims[1]/factor, dims[2]/factor), N))
	else
		pts = shape.points
	end
	angle = isa(shape.angle, Distribution) ? rand(shape.angle) : shape.angle
	return Blob(dims=dims, center=center, N=N, factor=factor, points=pts, angle=angle)
end

# ╔═╡ bc244e53-ddc4-4bf2-a6e4-74832e2e5e7d
rand(Ellipse(dims))

# ╔═╡ f8047ba9-727c-4f3f-9fda-04acef9f0683
etest = rand(Ellipse(dims=dims, center=[40,40], angle=float(π)))

# ╔═╡ 428b5b39-d83c-423b-811a-d660b0b53753
rand(etest)

# ╔═╡ 070e842b-5b9e-4631-adc4-0ec712d1008e
seed = 3

# ╔═╡ 00337f85-7b05-4a58-91ec-a65d257a992c
@bind d Slider([10, 20, 50, 80], default=80, show_value=true)

# ╔═╡ baa76cd0-4b68-4576-941a-e6f6fb0f7cf7
md"""
## Generate subsurface map
"""

# ╔═╡ f0ab6a97-20b8-4292-857c-c5059f6f7c81
function normalize01(mat)
	min = minimum(mat)
	max = maximum(mat)
	return (mat .- min) / (max - min)
end

# ╔═╡ 54407c87-7910-4225-bb12-2fc182ac3ed6
function Base.convert(::Type{Float64}, m::ColorTypes.ARGB32)
	return convert(Float64, Gray(m))
end

# ╔═╡ cce8da91-f857-4589-ac9b-640fd04a0934
function blur(img, σ)
	return imfilter(img, Kernel.gaussian(σ))
end

# ╔═╡ 299c3427-de7c-4380-9a2b-d496a63af9cc
function generate_subsurface(shape::Shape, σ=2)
	dims = shape.dims
	shape = rand(shape)
	mat = @imagematrix begin
		gsave()
		background("black")
		origin(Point(0,0))
		pos = Point(shape.center...)
		Luxor.translate(pos)
		if hasproperty(shape, :angle)
			Luxor.rotate(shape.angle)
		end
		grayscale = 0.4
		sethue(grayscale, grayscale, grayscale)
		draw(shape)
		grestore()
	end dims[1] dims[2]

	return convert(Matrix{Float64}, Gray.(blur(mat, σ)))
end

# ╔═╡ 496d4b3a-286d-4f6c-8c0f-4de0b7893ae1
begin
	Random.seed!(seed_ellipse)
	circle1 = rand(Circle(dims))
	circle2 = rand(Circle(dims))
	ore1_c_perfect = generate_subsurface(circle1, 0)
	ore2_c_perfect = generate_subsurface(circle2, 0)
	ore_c_perfect = ore1_c_perfect + ore2_c_perfect
	plot_subsurface(ore_c_perfect)
	title!("sampled shapes")
end

# ╔═╡ afe07e24-6bd3-4c5c-86bd-d746e22d4a44
begin
	ore1_c = generate_subsurface(circle1, 2)
	ore2_c = generate_subsurface(circle2, 2)
	ore_c = ore1_c + ore2_c
	plot_subsurface(ore_c)
	title!("blurred shapes")
end

# ╔═╡ 57720fcc-3272-4146-81a8-2f82c424814c
begin
	noise_c = imfilter(rand(dims...), Kernel.gaussian(1))
	plot_subsurface(noise_c, ore_c)
	title!("added background noise")
end

# ╔═╡ 3ae805bb-e5cf-4845-b407-96ea7c10fb31
begin
	ore_c_example =
		generate_subsurface(Circle(dims)) + generate_subsurface(Circle(dims))
	noise_c_example = imfilter(rand(dims...), Kernel.gaussian(1))
	plot_subsurface(noise_c_example, ore_c_example)
end

# ╔═╡ 8cf606e4-b897-472f-a8f4-3ab53a99f594
plot_subsurface(generate_subsurface(etest))

# ╔═╡ ca688c7e-f6ff-4b8c-93b0-d43aa0d4d940
begin
	Random.seed!(seed_ellipse)
	ellipse1 = rand(Ellipse(dims))
	ellipse2 = rand(Ellipse(dims))
	ore1_e_perfect = generate_subsurface(ellipse1, 0)
	ore2_e_perfect = generate_subsurface(ellipse2, 0)
	ore_e_perfect = ore1_e_perfect + ore2_e_perfect
	plot_subsurface(ore_e_perfect)
	title!("sampled shapes")
end

# ╔═╡ 0122ef4a-215e-4218-8b07-b0a3e8a42d96
begin
	ore1_e = generate_subsurface(ellipse1, 2)
	ore2_e = generate_subsurface(ellipse2, 2)
	ore_e = ore1_e + ore2_e
	plot_subsurface(ore_e)
	title!("blurred shapes")
end

# ╔═╡ 53311e80-c55d-46e3-ae49-7e7ab5808e79
begin
	noise_e = imfilter(rand(dims...), Kernel.gaussian(1))
	plot_subsurface(noise_e, ore_e)
	title!("added background noise")
end

# ╔═╡ 23180eef-c7ea-489f-b4f4-2dcc59796f5c
begin
	ore_e_example =
		generate_subsurface(Ellipse(dims)) + generate_subsurface(Ellipse(dims))
	noise_e_example = imfilter(rand(dims...), Kernel.gaussian(1))
	plot_subsurface(noise_e_example, ore_e_example)
end

# ╔═╡ da49bd52-857b-49f9-a7e4-09f24aa23b25
begin
	Random.seed!(seed)
	blob1 = rand(Blob(dims))
	blob2 = rand(Blob(dims))
	ore1_perfect = generate_subsurface(blob1, 0)
	ore2_perfect = generate_subsurface(blob2, 0)
	ore_perfect = ore1_perfect + ore2_perfect
	plot_subsurface(ore_perfect)
	title!("sampled shapes")
end

# ╔═╡ 1131d61e-bd4a-40d6-a07f-a60312ca52f2
begin
	ore1 = generate_subsurface(blob1, 2)
	ore2 = generate_subsurface(blob2, 2)
	ore = ore1 + ore2
	# noise = imfilter(rand(dims...), Kernel.gaussian(1))
	plot_subsurface(ore)
	title!("blurred shapes")
end

# ╔═╡ 8282e0e5-cd61-4526-a471-356b44e336ce
begin
	noise = imfilter(rand(dims...), Kernel.gaussian(1))
	plot_subsurface(noise, ore)
	title!("added background noise")
end

# ╔═╡ aaeb2772-16a9-4ec7-b8d6-37feb7c3adf2
begin
	ore_b_example =
		generate_subsurface(Blob(dims)) + generate_subsurface(Blob(dims))
	noise_b_example = imfilter(rand(dims...), Kernel.gaussian(1))
	plot_subsurface(noise_b_example, ore_b_example)
end

# ╔═╡ 1b66e642-f273-4bfa-9a72-a26eec713c38
begin
	scale_dims = (d, d)
	ground = imresize(noise_b_example + ore_b_example, scale_dims) |>
		p->plot_subsurface(zeros(scale_dims...), p)
	plot!(colorbar=false)
end

# ╔═╡ cce3cce1-6a04-4058-94a1-e9339a89ebe1
md"""
## Generate directly in low-fidelity
"""

# ╔═╡ 217f52eb-c926-4e25-86fe-053d6dccec3b
begin
	Random.seed!(0)
	dims_lf = (80, 80)
	shape1_lf = rand(Circle(dims_lf))
	shape2_lf = rand(Circle(dims_lf))
	ore1_e_lf = generate_subsurface(shape1_lf, 1)
	ore2_e_lf = generate_subsurface(shape2_lf, 1)
	ore_e_lf = ore1_e_lf + ore2_e_lf
	noise_lf = imfilter(rand(dims_lf...), Kernel.gaussian(exp(1/dims_lf[1])))
	plot_subsurface(noise_lf, ore_e_lf)
end

# ╔═╡ ccfeaa9c-4585-491f-b065-f979499bdd99
exp(1/80)

# ╔═╡ c1344ee7-7e0d-40e4-85da-26d70bc64ca1
md"""
## Particle filtering with rejection
"""

# ╔═╡ 46200976-234a-4d5c-8b0c-b38971c9dda8
pf_dims = (40,40)

# ╔═╡ 16461d4b-e8c2-41d6-8166-37e6f3acf89b
obs_pos = [9, 18]

# ╔═╡ dce1bd6d-4e13-460b-a7fb-8478c5aed8ec
obs_pos2 = [12, 19]

# ╔═╡ 001c7316-a0b0-45b2-b25d-ee14b6e86cfe
obs_pos3 = [22, 16] # [11, 16]

# ╔═╡ 240e4a06-ca1b-46fb-b3e8-bf59eb913e66
obs_pos4 = [25, 15] # [21, 9] # [11, 21] # [18, 12]

# ╔═╡ 1bb50375-3a64-4651-91c4-d790cdf41fec
@bind temp_seed Slider(1:100, show_value=true, default=29)

# ╔═╡ 4b8ff8f4-3c44-43b8-a463-d63a9126189e
# plot([begin
# 	plot_subsurface(G(particles))
# 	plot!(colorbar=false, axis=[])
# end for _ in 1:9]...)

# ╔═╡ 4900d76a-b512-49c6-9bd4-d26c52d9f5cc
function rejection_filtering(G, obs_pos, obs_val; m=20, δ=0.05, max_samples=100_000)
	particles = []
	weights = []
	for _ in 1:max_samples
		length(particles) ≥ m && break
		sample = G() # generate new sample
		sample_val = sample[obs_pos...]
		# accept that sample if it's consistent (within some threshold δ)
		difference = abs(obs_val - sample_val)
		if difference ≤ δ
			push!(particles, sample)
			push!(weights, sample_val/difference) # NOTE: 1/difference
		end
	end
	return particles, normalize(weights, 1)
end

# ╔═╡ c99caffd-1cad-43b9-a087-8944c31e597f
function conditional_rejection_filtering(G, particles, obs_pos, obs_val,
		prev_obs_pos, prev_obs_val; m=20, δ=0.05, max_samples=10_000)
	conditioned_particles = []
	for sample in particles
		length(conditioned_particles) ≥ m && break
		sample_val = sample[obs_pos...]
		if abs(obs_val - sample_val) ≤ δ
			push!(conditioned_particles, sample)
		end
	end
	if length(conditioned_particles) < m
		@warn "Not enough"
		# resample until we get m conditioned particles
		# for _ in 1:max_samples
			
		# for _ in 1:max_samples
		# 	length(conditioned_particles) ≥ m && break
		# 	sample = G(particles) # generate new sample, seeding with sampled center
		# 	sample_val = sample[obs_pos...]
		# 	if abs(obs_val - sample_val) ≤ δ
		# 		@info "–"
		# 		push!(conditioned_particles, sample)
		# 	end
		# end
	end
	return conditioned_particles
end

# ╔═╡ 910cae01-75f6-47d7-803a-8d54468821e5
# function conditional_rejection_filtering_core()

# ╔═╡ 8256e288-708f-451f-88d8-af8f29d136e6
function plot_samples(particles, obs_pos=nothing)
	particle_plots = [begin
		plot_subsurface(p)
		plot!(colorbar=false, axis=[])
		if !isnothing(obs_pos)
			scatter!([obs_pos[2]], [obs_pos[1]], c=:red, label=false, ms=2)
		else
			plot!()
		end		
	end for p in particles]
	plot(particle_plots...)
end

# ╔═╡ aae4a69e-6bdf-4c98-9a40-733a89e47625
function plot_volume(particles)
	volumes = [sum(p .≥ massive_threshold) for p in particles]
	μ = round(mean(volumes), digits=2)
	σ = round(std(volumes), digits=2)
	h = fit(Distributions.Histogram, volumes, [0:10:300;])
	h = normalize(h, mode=:probability)
	height = maximum(h.weights)
	plot(h, label=false, c=:gray, ylims=(0, height*1.05))
	plot!([μ, μ], [0, height], c=:red, lw=3, label="mean")
	title!("belief volumes: μ = $μ, σ = $σ")	
end

# ╔═╡ 6e01aa35-784d-483b-9404-fbc1caf9244a
function weighted_mean(values, weights)
	return sum(weights[i] * values[i] for i in 1:length(values)) / sum(weights)
end

# ╔═╡ 8a6acb2a-d7bc-4d98-9e7a-c110888ef0da
@bind is_weighted CheckBox()

# ╔═╡ cb258b50-613d-4ac7-b6ae-8b795b318f60
md"""
## MCMC
- Use candidate "centers" to resample starting from that center and continue to fill the `conditioned_particle` set (basing it off of the `conditioned_particles` because we know those are in agreement)
"""

# ╔═╡ b0fea9ee-17d8-4e9d-90bc-b7a7d916430f
domain = CartesianIndices(pf_dims)

# ╔═╡ 8aa3ad9d-b4e2-4452-bef3-fb4d8916ed77
function round2dims(x, dims)
	x = round.(Int, x)
	x[1] = clamp(x[1], 1, dims[1])
	x[2] = clamp(x[2], 1, dims[2])
	return x
end

# ╔═╡ 9d2cc2be-b162-41da-99a9-501f47d414bf
function mcmc(f, startx, domain, g=x->MvNormal(x,[0.01 0; 0 0.01]);
		T=200, seed=nothing)
	!isnothing(seed) && Random.seed!(seed)
	x = x′ = startx
	if isa(domain, CartesianIndices)
		x = [x.I...]
		x′ = [x′.I...]
	end
	y = f(x)
	cost = -Inf
	X::Vector{Vector{Real}} = [x]
	Y::Vector{Real} = [y]
	for t in 1:T
		x = rand(g(x′)) # sample point close to previous accepted point
		x = round2dims(x, size(domain)) # restrict to grid
		while x ∉ X
			x = rand(g(x′))
			x = round2dims(x, size(domain))
		end
		y = f(x) # evaluate point
		Q = pdf(g(x′), x) # g(x | x′)
		Q′ = pdf(g(x), x′) # g(x′ | x)
		α = min(1, y/cost * Q′/Q)
		if α == 1
			x′, cost = x, y # accept
			push!(X, x); push!(Y, y)
		else
			cost = rand() > α ? y : cost # accept with α probability ∝ difference
		end
	end

	return X, Y, X[argmax(Y)]
end

# ╔═╡ b569794a-0bc0-4f04-b366-1593878aae5d
md"""
```julia
for t in 1:T
	x = rand(g(x′)) # sample point close to previous accepted point
	y = f(x) # evaluate point
	if y > cost
		x′, cost = x, y # accept
		push!(X, x); push!(Y, y)
	else
		α = y/cost
		cost = rand() > α ? y : cost # accept with α probability ∝ difference
	end
end
```"""

# ╔═╡ 5fee2e36-b409-4b42-9c95-3e1eea65eb1b
function plot_mcmc(xy; trace=true)
	if isa(xy, Vector)
		x = first.(xy)
		y = last.(xy)
		c = :red # [get(grad, i/length(x)) for i in 1:length(x)]
		if trace
			# Notice: y, x flipped (row, col)
			plot!(y, x, label=false, lw=1/2, ms=2, m=:circle, c=c)
		else
			# Notice: y, x flipped (row, col)
			scatter!(y, x, c=:red, label=false, ms=2)
		end
	end
end

# ╔═╡ 81b3e18d-2627-4ee2-945c-e4a92259e6b3
# plot(cat)

# ╔═╡ 3a937351-78ba-4a2f-b804-f27e96bdaebe
function normalize_particles(particles)
	A = vec(particles)
	normed_zero_one = A .- minimum(A)
	normed_sum_one = normalize(normed_zero_one, 1)
	return normed_sum_one
end

# ╔═╡ a2ad5bb1-c088-4da3-978e-b5c88b9c6463
function G(particles=nothing; shape=nothing, σ=1, with_noise=true)
	if isnothing(particles)
		if isnothing(shape)
			shape = Ellipse(pf_dims)
		end
	else
		μ_particles = mean(particles)
		# turn mean field into PMF
		distr = Categorical(normalize_particles(μ_particles))
		center_sample = rand(distr) # sample for a new center (heuristic)
		center_coord = [CartesianIndices(size(μ_particles))[center_sample].I...]
		shape = Ellipse(dims=pf_dims, center=center_coord)
	end
	ore = generate_subsurface(shape, σ)
	if with_noise
		noise = imfilter(rand(pf_dims...), Kernel.gaussian(2))
		return ore + noise
	else
		return ore
	end
end

# ╔═╡ f974450e-800a-4547-b4a8-935dd77bc54f
if false
	Random.seed!(temp_seed)
	temp_truth = G(shape=Blob(pf_dims), σ=0, with_noise=false)
	plot_subsurface(temp_truth)
	title!("truth without noise/blur")
end

# ╔═╡ dc2f44ed-248f-4c8a-a124-22c73c59532e
begin
	Random.seed!(29)
	truth = G(shape=Blob(pf_dims))
	plot_subsurface(truth)
	scatter!([obs_pos[2]], [obs_pos[1]], c=:red, label=false)
	scatter!([obs_pos2[2]], [obs_pos2[1]], c=:orange, label=false)
	scatter!([obs_pos3[2]], [obs_pos3[1]], c=:yellow, label=false)
	scatter!([obs_pos4[2]], [obs_pos4[1]], c=:white, label=false)
	title!("truth w/ observations")
end

# ╔═╡ 3c02560f-dc60-4c38-a765-46bfd204a25b
obs_val = truth[obs_pos...]

# ╔═╡ 2dff5f86-cc27-4da7-b326-3af959714949
obs_val/0.05

# ╔═╡ faa23422-2246-4c95-8fd3-b9e1658c546f
obs_val/0.005

# ╔═╡ ae024dfa-4b3f-42dd-b43d-e4654345ebdb
obs_val2 = truth[obs_pos2...]

# ╔═╡ 7dcc0530-89f8-4780-a812-6a085db2efb9
obs_val3 = truth[obs_pos3...]

# ╔═╡ 1e64a718-bf1f-4113-80a7-fda421491f00
obs_val4 = truth[obs_pos4...]

# ╔═╡ 7305f3e9-d8ba-4bc5-a6a9-817878f965d3
norm_truth = (truth .- minimum(truth)) / (maximum(truth) - minimum(truth));

# ╔═╡ 74207c6f-7b17-451e-8127-f4ac99f4c751
plot(
	plot_subsurface(truth .≥ 0.7),
	plot_subsurface(norm_truth .≥ 0.7)
)

# ╔═╡ 7ed162db-379a-460c-a6e8-04021e267d46
sum(truth .≥ massive_threshold), sum(norm_truth .≥ massive_threshold)

# ╔═╡ 0b378769-ac0b-4d62-a8c0-d19cf1bef672
minimum(norm_truth), maximum(norm_truth)

# ╔═╡ db8d8b36-8208-4a25-9240-92fd75afefa1
truth |> plot_subsurface

# ╔═╡ 3a39721e-323b-4726-8ca8-de894a6103ab
particles, weights = rejection_filtering(G, obs_pos, obs_val; m=1000, δ=0.05);

# ╔═╡ 76e8280f-8899-4ed9-8e68-a215ec1a9685
length(particles)

# ╔═╡ 4a209db9-53f9-4f06-8973-78430009badc
plot_samples(particles[1:6], obs_pos)

# ╔═╡ 339f8dad-728d-456c-ab6d-8e210b42e1e4
plot(begin
		plot_subsurface(mean(particles), clims=false)
		scatter!([obs_pos[2]], [obs_pos[1]], c=:red, label=false)
		title!("belief mean")
	end,
	begin
		plot_subsurface(std(particles), clims=false)
		title!("belief std")
	end)

# ╔═╡ 8df8f92e-859e-432c-ace0-552ffae42540
plot_volume(particles)

# ╔═╡ daff6450-e382-4f5b-b6bb-152899efef4c
plot(weights)

# ╔═╡ 3d2328b7-b746-413c-8b70-8a0e84ef620e
if is_weighted
	weighted_mean(particles, weights) |> x->plot_subsurface(x, clims=false)
else
	mean(particles) |> x->plot_subsurface(x, clims=false)
	# truth |> plot_subsurface
end

# ╔═╡ 0758cf10-2e61-446f-8913-a99f9401cb13
vec(mean(particles))[123]

# ╔═╡ 96d75c69-e56e-4084-9b47-f38bd784dd9e
mean(particles)[123]

# ╔═╡ 3b8f91ff-6dcd-44c1-8aee-897c7a0cfc1c
f = x->begin
	μ_particles = weighted_mean(particles, weights)
	# μ_particles = mean(particles)
	if isa(x, CartesianIndex)
		coord = x.I
	else
		coord = round2dims(x, size(μ_particles))
	end
	return μ_particles[coord[1], coord[2]]
end

# ╔═╡ f80b6d79-7930-4a1b-9abc-4f6e86e5eaf8
var_var = 1e-4/var(var(particles))

# ╔═╡ 3c58bb53-2fae-4ba9-b331-e204438950da
g = x->MvNormal(x, [var_var 0; 0 var_var])

# ╔═╡ deb8157f-1c7f-4326-9571-30dc12396826
startx = argmax(weighted_mean(particles, weights))

# ╔═╡ 9f0584cb-3887-4a6b-9905-c1085fc5c736
X_mcmc, Y_mcmc, Xarg_mcmc = mcmc(f, startx, domain, g; T=1000)

# ╔═╡ e820625d-e591-4da8-bced-5765c88e7b59
length(X_mcmc), length(unique(X_mcmc))

# ╔═╡ d4fb8df3-62ce-4216-80ed-d14d2e52662e
begin
	plot_subsurface(weighted_mean(particles, weights), clims=false)
	# plot_subsurface(mean(particles), clims=false)
	plot_mcmc(X_mcmc)
end

# ╔═╡ 0854cb6c-3f27-45a0-ae2e-bf1e8fc647f8
cat = Categorical(normalize(vec(mean(particles)), 1))

# ╔═╡ 39425255-756c-475a-81a5-0fbf334847c0
begin
	plot_subsurface(normalize(weighted_mean(particles, weights),1), clims=false)
	for _ in 1:100
		# distr = Categorical(normalize_particles(mean(particles)))
		distr = Categorical(normalize_particles(weighted_mean(particles, weights)))
		center_sample = rand(distr)		
		center_coord = CartesianIndices(size(mean(particles)))[center_sample]
		scatter!([center_coord[2]], [center_coord[1]], c=:red, label=false, alpha=0.5)
	end
	plot!()
end

# ╔═╡ 63653c20-13b3-4fab-98f3-b1a8aff80d83
cat_normed = Categorical(normalize_particles(mean(particles)))

# ╔═╡ 3ce20506-e600-4641-9ea6-5bff5bb2c94e
# plot(cat_normed)

# ╔═╡ 527ebc79-8a8c-4a72-a7e1-1f584078a88a
plot_subsurface(reshape(normalize_particles(mean(particles)), pf_dims), clims=false)

# ╔═╡ 0ca7cd5c-cd98-4c16-8fa9-e8b7caf35850
md"""
## Conditioning
- turn mean field into probability distribution
    - sample center point from mean field
"""

# ╔═╡ b1545c2e-53f3-4e79-afb4-c5b3b44a67c4
cond_particles =
	conditional_rejection_filtering(G, particles, obs_pos2, obs_val2, obs_pos, obs_val; m=1000);

# ╔═╡ 3b8439d6-5295-4ce3-9d78-a317a8fe76b2
length(cond_particles)

# ╔═╡ a68987db-37a9-474d-a98d-89833d926801
plot_samples(cond_particles[1:6], obs_pos2)

# ╔═╡ d6a33876-f825-4e69-be15-33bd670ab152
plot(begin
		plot_subsurface(mean(cond_particles), clims=false)
		scatter!([obs_pos2[2]], [obs_pos2[1]], c=:orange, label=false)
	end,
	plot_subsurface(std(cond_particles), clims=false))

# ╔═╡ 010d5830-2949-4b44-9fda-b2094a7cfdfc
plot_volume(cond_particles)

# ╔═╡ 384d8315-6f2a-4c3f-aa22-d23260aece37
md"""
### Conditioning (2)
"""

# ╔═╡ 2c751941-f20f-4844-aea6-8f7553ebd525
cond_particles2 =
	conditional_rejection_filtering(G, cond_particles, obs_pos3, obs_val3, obs_pos2, obs_val2; m=1000);

# ╔═╡ 9afeccaf-85b2-457a-996f-086b8684568f
length(cond_particles2)

# ╔═╡ 67e96163-93e9-4882-8948-730a386b3208
plot_samples(cond_particles2[1:min(length(cond_particles), 6)], obs_pos3)

# ╔═╡ dba0ec68-701e-4c60-9d0a-2bd44563ac7c
plot(begin
		plot_subsurface(mean(cond_particles2), clims=false)
		scatter!([obs_pos3[2]], [obs_pos3[1]], c=:yellow, label=false)
	end,
	plot_subsurface(std(cond_particles2), clims=false))

# ╔═╡ 3a9b4d7f-4387-463e-af33-d090d735c015
plot_volume(cond_particles2)

# ╔═╡ 96aa6311-a384-422a-8b14-3f186c2c0e25
md"""
### Conditioning (3)
"""

# ╔═╡ 874a4572-cdd1-44f5-a4ff-71fcd8079914
cond_particles3 =
	conditional_rejection_filtering(G, cond_particles2, obs_pos4, obs_val4, obs_pos3, obs_val3; m=1000);

# ╔═╡ 97c065b4-6d41-49fc-9db3-ac553d08db83
length.([particles,cond_particles, cond_particles2, cond_particles3])

# ╔═╡ c3702315-3d32-4a01-89bd-299b3d278c21
plot_subsurface(mean(cond_particles3) .≥ massive_threshold)

# ╔═╡ c00aa719-e73f-4f33-a72f-c85750056d05
length(cond_particles3)

# ╔═╡ 1b9fe5c2-a96f-4d68-9738-deb774627fa1
plot_samples(cond_particles3[1:min(length(cond_particles3), 6)], obs_pos4)

# ╔═╡ 630353a1-93cc-47e5-80ac-7c6cf5b10583
plot(begin
		plot_subsurface(mean(cond_particles3), clims=false)
		scatter!([obs_pos4[2]], [obs_pos4[1]], c=:white, label=false)
		title!("mean")
	end,
	begin
		plot_subsurface(std(cond_particles3), clims=false)
		title!("std")
	end)

# ╔═╡ 57c2fc19-b968-4e99-b866-5605507da2c3
plot_volume(cond_particles3)

# ╔═╡ 085ad6cb-022a-4eb7-8a6f-06c6f93fd48b
argmax(mean(cond_particles3))

# ╔═╡ bfef858b-c24f-4839-b8f4-31ac5f014401
plot(
	plot_subsurface(truth .≥ massive_threshold,
		title=string("(truth) massive ore = ", sum(truth .≥ massive_threshold))),
	plot_subsurface(mean(cond_particles3) .≥ massive_threshold,
		title=string("(belief) massive ore = ", sum(mean(cond_particles3) .≥ massive_threshold)))
)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
ImageFiltering = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
POMDPs = "a93abf59-7444-517b-a68a-c42f96afdd7d"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"

[compat]
Distributions = "~0.25.41"
ImageFiltering = "~0.7.1"
Images = "~0.25.1"
Luxor = "~3.0.0"
POMDPs = "~0.9.3"
Parameters = "~0.12.3"
Plots = "~1.25.7"
PlutoUI = "~0.7.32"
ProgressLogging = "~0.1.4"
Revise = "~3.3.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
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

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "ffc6588e17bcfcaa79dfa5b4f417025e755f83fc"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "d127d5e4d86c7680b20c35d40b503c74b9a39b5e"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.4"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

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

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "9aa8a5ebb6b5bf469a7e0e2b5202cf6f8c291104"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.6"

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

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

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

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

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

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "5863b0b10512ed4add2b5ec07e335dc6121065a5"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.41"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "84f04fe68a3176a583b864e492578b9466d87f1e"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

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

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "67551df041955cc6ee2ed098718c8fcd7fc7aebe"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.12.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

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
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

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

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "d727758173afef0af878b29ac364a0eca299fc6b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.5.1"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageContrastAdjustment]]
deps = ["ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "0d75cafa80cf22026cea21a8e6cf965295003edc"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.10"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "7a20463713d239a19cbad3f6991e404aca876bda"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.15"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "15bd05c1c0d5dbb32a9a3d7e0ad2d50dd6167189"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.1"

[[deps.ImageIO]]
deps = ["FileIO", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "816fc866edd8307a6e79a575e6585bfab8cef27f"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.0"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f025b79883f361fa1bd80ad132773161d231fd9f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.12+2"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.ImageMorphology]]
deps = ["ImageCore", "LinearAlgebra", "Requires", "TiledIteration"]
git-tree-sha1 = "7668b123ecfd39a6ae3fc31c532b588999bdc166"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.3.1"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "OffsetArrays", "Statistics"]
git-tree-sha1 = "1d2d73b14198d10f7f12bf7f8481fd4b3ff5cd61"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.0"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "36832067ea220818d105d718527d6ed02385bf22"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.7.0"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "d0ac64c9bee0aed6fdbb2bc0e5dfa9a3a78e3acc"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.3"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "42fe8de1fe1f80dab37a39d391b6301f7aeaa7b8"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.9.4"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "11d268adba1869067620659e7cdf07f5e54b6c76"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.25.1"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "00019244715621f473d399e4e1842e479a69a42e"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

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

[[deps.JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "bcb31db46795eeb64480c89d854615bc78a13289"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.19"

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
git-tree-sha1 = "6ca01d8e5bc75d178e8ac2d1f741d02946dc1853"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.2"

[[deps.Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

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

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "25d5e6b4eb3558613ace1c67d6a871420bfca527"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.52.4+0"

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

[[deps.LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

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

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "Dates", "FFMPEG", "FileIO", "Juno", "LaTeXStrings", "Random", "Requires", "Rsvg"]
git-tree-sha1 = "81a4fd2c618ba952feec85e4236f36c7a5660393"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "3.0.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

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

[[deps.Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "2af69ff3c024d13bde52b34a2a7d6887d4e7b438"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[deps.NamedTupleTools]]
git-tree-sha1 = "63831dcea5e11db1c0925efe5ef5fc01d528c522"
uuid = "d9ec5142-1e00-5aa0-9d6a-321866360f50"
version = "0.13.7"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "16baacfdc8758bc374882566c9187e785e85c2f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.9"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

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

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "6d105d40e30b635cfed9d52ec29cf456e27d38f8"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.12"

[[deps.POMDPLinter]]
deps = ["Logging"]
git-tree-sha1 = "cee5817d06f5e1a9054f3e1bbb50cbabae4cd5a5"
uuid = "f3bd98c0-eb40-45e2-9eb1-f2763262d755"
version = "0.1.1"

[[deps.POMDPs]]
deps = ["Distributions", "LightGraphs", "NamedTupleTools", "POMDPLinter", "Pkg", "Random", "Statistics"]
git-tree-sha1 = "3a8f6cf6a3b7b499ec4294f2eb2b16b9dc8a7513"
uuid = "a93abf59-7444-517b-a68a-c42f96afdd7d"
version = "0.9.3"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9bc1871464b12ed19297fbc56c4fb4ba84988b0d"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.47.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

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
git-tree-sha1 = "ae6145ca68947569058866e443df69587acc1806"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.32"

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

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

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

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra"]
git-tree-sha1 = "adf644ef95a5e26c8774890a509a55b7791a139f"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "01d341f502250e81f6fec0afe662aa861392a3aa"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.2"

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

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

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

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "2f9d4d6679b5f0394c52731db3794166f49d5131"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.3.1"

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

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "405148000e80f70b31e7732ea93288aecb1793fa"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.2.0"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "3d3dc66eb46568fb3a5259034bfc752a0eb0c686"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.0.0"

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

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "a6f404cc44d3d3b28c793ec0eb59af709d827e4e"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.2.1"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

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

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "b4912cd034cdf968e06ca5f943bb54b17b97793a"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2884859916598f974858ff01df7dfc6c708dd895"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.3"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "991d34bbff0d9125d93ba15887d6594e8e84b305"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.3"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

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

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

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

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "c23323cd30d60941f8c68419a70905d9bdd92808"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.6+1"

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

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

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
# ╟─668a9810-7a3f-11ec-0861-89fc669dae9e
# ╠═909a415e-f420-441c-bd2d-4ec86239113b
# ╠═7fb2c2da-ca38-489c-b8c1-1cc402b68222
# ╠═7062f835-165c-414b-b25f-a060de10b3cc
# ╠═58dab9b5-25b2-42fe-b918-0dea561d5fd5
# ╠═a4aa5758-6caf-47eb-a18c-ff2281453e23
# ╠═4d70219e-9fd3-4d14-80d2-3c12fff9e837
# ╠═bd1330a9-b65a-46a5-93ce-ecaf29632699
# ╟─809bb5c7-8635-451a-8bba-11c3b5da215b
# ╠═b216e135-e2a2-4646-9df4-1ccdda7d2dec
# ╠═1ed349dc-a7d1-40b9-8c2b-e915d72fd150
# ╠═a91efc56-40e2-4a40-af75-b60d06f933c6
# ╠═28dd3755-9c37-4e30-b58f-4b04c3ded69e
# ╠═cee46dc3-32d4-4e79-96c4-5480f1402b80
# ╟─929486d5-1e98-4e03-a2d6-3279699d95df
# ╠═718b84cd-9851-4b4a-b9fb-b0d1ec2a628e
# ╠═e78b494c-d32c-4123-a19c-9b1ec2fca78b
# ╠═e26ba54c-efd3-44af-b8df-79569a31c160
# ╠═ee6b7a25-2232-4130-81ab-5ddcf94e0b4f
# ╠═14e5bbf4-f51b-4d7a-9eba-cbe71b99fe61
# ╠═355aebe8-a750-434e-8344-1d8b5ecc33fe
# ╠═8113c944-287f-442f-8c43-a00fb39b9ca5
# ╟─a7d933d7-2e99-4303-bf8a-7b57981f5349
# ╠═ce783848-56cc-43d6-907c-74dfd6d8c0cd
# ╠═f395c41f-399f-48f6-9b8c-3010d555887d
# ╠═3d8da47d-bee1-46e8-96ed-07e1c5ad897d
# ╠═8f5e6e85-83a8-4410-8f8f-d84b6925b549
# ╠═496d4b3a-286d-4f6c-8c0f-4de0b7893ae1
# ╠═afe07e24-6bd3-4c5c-86bd-d746e22d4a44
# ╠═57720fcc-3272-4146-81a8-2f82c424814c
# ╠═3ae805bb-e5cf-4845-b407-96ea7c10fb31
# ╟─3dd84717-7214-4d70-b82d-d1e0558c129a
# ╠═bc48337a-8699-45a5-a135-de75933e07a4
# ╠═327409bc-bbf4-4236-9de5-19a61b69b9a1
# ╠═06e993f6-4300-479b-8f4f-3c1666f29acc
# ╠═bc244e53-ddc4-4bf2-a6e4-74832e2e5e7d
# ╠═f8047ba9-727c-4f3f-9fda-04acef9f0683
# ╠═428b5b39-d83c-423b-811a-d660b0b53753
# ╠═8cf606e4-b897-472f-a8f4-3ab53a99f594
# ╠═63799459-b343-4162-a421-ba3a77d50b65
# ╠═ca688c7e-f6ff-4b8c-93b0-d43aa0d4d940
# ╠═0122ef4a-215e-4218-8b07-b0a3e8a42d96
# ╠═53311e80-c55d-46e3-ae49-7e7ab5808e79
# ╠═23180eef-c7ea-489f-b4f4-2dcc59796f5c
# ╟─8fc9b648-3810-4b32-8da8-7bf14225c6e8
# ╠═2014a964-ab02-4fc4-962f-2f057a1da115
# ╠═9e7ca98c-34a2-44bb-ad72-e66e12dc006d
# ╠═336f0ea3-ea05-4207-bf2a-e21ac380f424
# ╠═070e842b-5b9e-4631-adc4-0ec712d1008e
# ╠═da49bd52-857b-49f9-a7e4-09f24aa23b25
# ╠═1131d61e-bd4a-40d6-a07f-a60312ca52f2
# ╠═8282e0e5-cd61-4526-a471-356b44e336ce
# ╠═aaeb2772-16a9-4ec7-b8d6-37feb7c3adf2
# ╠═e487916b-3581-4814-b2cd-36939fb471c7
# ╠═00337f85-7b05-4a58-91ec-a65d257a992c
# ╠═1b66e642-f273-4bfa-9a72-a26eec713c38
# ╟─baa76cd0-4b68-4576-941a-e6f6fb0f7cf7
# ╠═299c3427-de7c-4380-9a2b-d496a63af9cc
# ╠═f0ab6a97-20b8-4292-857c-c5059f6f7c81
# ╠═54407c87-7910-4225-bb12-2fc182ac3ed6
# ╠═cce8da91-f857-4589-ac9b-640fd04a0934
# ╟─cce3cce1-6a04-4058-94a1-e9339a89ebe1
# ╠═217f52eb-c926-4e25-86fe-053d6dccec3b
# ╠═ccfeaa9c-4585-491f-b065-f979499bdd99
# ╟─c1344ee7-7e0d-40e4-85da-26d70bc64ca1
# ╠═522efd58-3756-4ad7-931e-ef1f0bd8c483
# ╠═46200976-234a-4d5c-8b0c-b38971c9dda8
# ╠═16461d4b-e8c2-41d6-8166-37e6f3acf89b
# ╠═3c02560f-dc60-4c38-a765-46bfd204a25b
# ╠═dce1bd6d-4e13-460b-a7fb-8478c5aed8ec
# ╠═ae024dfa-4b3f-42dd-b43d-e4654345ebdb
# ╠═001c7316-a0b0-45b2-b25d-ee14b6e86cfe
# ╠═7dcc0530-89f8-4780-a812-6a085db2efb9
# ╠═240e4a06-ca1b-46fb-b3e8-bf59eb913e66
# ╠═1e64a718-bf1f-4113-80a7-fda421491f00
# ╠═1bb50375-3a64-4651-91c4-d790cdf41fec
# ╠═f974450e-800a-4547-b4a8-935dd77bc54f
# ╠═dc2f44ed-248f-4c8a-a124-22c73c59532e
# ╠═97c065b4-6d41-49fc-9db3-ac553d08db83
# ╠═74207c6f-7b17-451e-8127-f4ac99f4c751
# ╠═7ed162db-379a-460c-a6e8-04021e267d46
# ╠═7305f3e9-d8ba-4bc5-a6a9-817878f965d3
# ╠═0b378769-ac0b-4d62-a8c0-d19cf1bef672
# ╠═a2ad5bb1-c088-4da3-978e-b5c88b9c6463
# ╠═4b8ff8f4-3c44-43b8-a463-d63a9126189e
# ╠═4900d76a-b512-49c6-9bd4-d26c52d9f5cc
# ╠═2dff5f86-cc27-4da7-b326-3af959714949
# ╠═faa23422-2246-4c95-8fd3-b9e1658c546f
# ╠═c99caffd-1cad-43b9-a087-8944c31e597f
# ╠═910cae01-75f6-47d7-803a-8d54468821e5
# ╠═3a39721e-323b-4726-8ca8-de894a6103ab
# ╠═76e8280f-8899-4ed9-8e68-a215ec1a9685
# ╠═8256e288-708f-451f-88d8-af8f29d136e6
# ╠═4a209db9-53f9-4f06-8973-78430009badc
# ╠═339f8dad-728d-456c-ab6d-8e210b42e1e4
# ╠═8df8f92e-859e-432c-ace0-552ffae42540
# ╠═aae4a69e-6bdf-4c98-9a40-733a89e47625
# ╠═c3702315-3d32-4a01-89bd-299b3d278c21
# ╠═daff6450-e382-4f5b-b6bb-152899efef4c
# ╠═6e01aa35-784d-483b-9404-fbc1caf9244a
# ╠═8a6acb2a-d7bc-4d98-9e7a-c110888ef0da
# ╠═3d2328b7-b746-413c-8b70-8a0e84ef620e
# ╠═db8d8b36-8208-4a25-9240-92fd75afefa1
# ╠═0758cf10-2e61-446f-8913-a99f9401cb13
# ╠═96d75c69-e56e-4084-9b47-f38bd784dd9e
# ╠═2039b5c2-b834-495a-a3dd-ea20394f9bef
# ╠═cb258b50-613d-4ac7-b6ae-8b795b318f60
# ╠═b0fea9ee-17d8-4e9d-90bc-b7a7d916430f
# ╠═8aa3ad9d-b4e2-4452-bef3-fb4d8916ed77
# ╠═3b8f91ff-6dcd-44c1-8aee-897c7a0cfc1c
# ╠═f80b6d79-7930-4a1b-9abc-4f6e86e5eaf8
# ╠═3c58bb53-2fae-4ba9-b331-e204438950da
# ╠═9d2cc2be-b162-41da-99a9-501f47d414bf
# ╟─b569794a-0bc0-4f04-b366-1593878aae5d
# ╠═deb8157f-1c7f-4326-9571-30dc12396826
# ╠═9f0584cb-3887-4a6b-9905-c1085fc5c736
# ╠═e820625d-e591-4da8-bced-5765c88e7b59
# ╠═d4fb8df3-62ce-4216-80ed-d14d2e52662e
# ╠═5fee2e36-b409-4b42-9c95-3e1eea65eb1b
# ╠═39425255-756c-475a-81a5-0fbf334847c0
# ╠═0854cb6c-3f27-45a0-ae2e-bf1e8fc647f8
# ╠═81b3e18d-2627-4ee2-945c-e4a92259e6b3
# ╠═3a937351-78ba-4a2f-b804-f27e96bdaebe
# ╠═63653c20-13b3-4fab-98f3-b1a8aff80d83
# ╠═3ce20506-e600-4641-9ea6-5bff5bb2c94e
# ╠═527ebc79-8a8c-4a72-a7e1-1f584078a88a
# ╠═0ca7cd5c-cd98-4c16-8fa9-e8b7caf35850
# ╠═b1545c2e-53f3-4e79-afb4-c5b3b44a67c4
# ╠═3b8439d6-5295-4ce3-9d78-a317a8fe76b2
# ╠═a68987db-37a9-474d-a98d-89833d926801
# ╠═d6a33876-f825-4e69-be15-33bd670ab152
# ╠═010d5830-2949-4b44-9fda-b2094a7cfdfc
# ╟─384d8315-6f2a-4c3f-aa22-d23260aece37
# ╠═2c751941-f20f-4844-aea6-8f7553ebd525
# ╠═9afeccaf-85b2-457a-996f-086b8684568f
# ╠═67e96163-93e9-4882-8948-730a386b3208
# ╠═dba0ec68-701e-4c60-9d0a-2bd44563ac7c
# ╠═3a9b4d7f-4387-463e-af33-d090d735c015
# ╟─96aa6311-a384-422a-8b14-3f186c2c0e25
# ╠═874a4572-cdd1-44f5-a4ff-71fcd8079914
# ╠═c00aa719-e73f-4f33-a72f-c85750056d05
# ╠═1b9fe5c2-a96f-4d68-9738-deb774627fa1
# ╠═630353a1-93cc-47e5-80ac-7c6cf5b10583
# ╠═57c2fc19-b968-4e99-b866-5605507da2c3
# ╠═085ad6cb-022a-4eb7-8a6f-06c6f93fd48b
# ╠═bfef858b-c24f-4839-b8f4-31ac5f014401
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
