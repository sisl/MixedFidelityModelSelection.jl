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

# ╔═╡ c2cfe5bc-ac18-4310-8a30-6b89c0a1e9d8
begin
	using Revise
	using Random
end

# ╔═╡ 9408fe18-99c9-48ff-9f6a-63caa2c99d96
using Plots; default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

# ╔═╡ 539de988-9349-47b1-8129-d54d80663154
using POMDPs

# ╔═╡ e8b732da-65b3-4977-879d-87f256ec48c1
using Parameters

# ╔═╡ 0e1a8833-0ebd-4893-a463-7dc092a60540
using Distributions

# ╔═╡ 6aba071e-5b44-43ae-b94c-1ad5bcd8baba
using Luxor

# ╔═╡ e5bc2147-fcdf-449e-9693-3d902e27bbc1
using ImageFiltering

# ╔═╡ 2baf8598-7964-4e77-9a6f-94427cf347bc
using PlutoUI

# ╔═╡ fdef4b60-7fbc-11ec-3adc-a7773dfccfa2
md"""
# Mainbody shapes
"""

# ╔═╡ fc4e7fc0-9b2c-4e3c-80c9-7120ef1019c5
massive_threshold = 0.7

# ╔═╡ 8a71c40c-2319-479e-8fb9-710a57d540b3
md"""
## Shapes
Non-circular shapes includes a rotation angle $\alpha$
- Square: side length $\ell$
- Circle: radius $r$
- Ellipse: width $w$ and height $h$
"""

# ╔═╡ feb16f09-aba5-444b-a814-1c6a816a6c2f
grid_dims = (50,50,1)

# ╔═╡ 37fc09d9-cd50-4197-879a-045b6af2e707
function limits(grid_dims)
	lims_x = (0, grid_dims[1]) .+ 1/2
	lims_y = (0, grid_dims[2]) .+ 1/2
	return (lims_x, lims_y)
end

# ╔═╡ 052eaac9-bac5-439a-8a56-058d81d7bb97
lims_x, lims_y = limits(grid_dims)

# ╔═╡ 14fdb72c-28b4-4aab-b55b-32ed943b8ed8
function plot_subsurface(noise, ore; clims=true, title="")
	grid_dims = size(ore)
	lims_x, lims_y = limits(grid_dims)
	subsurface = noise + ore
	heatmap(subsurface,
		aspect_ratio=1,
		xlims=lims_x,
		ylims=lims_y,
		c=:viridis,
		clims=clims ? (0,1) : (minimum(subsurface), maximum(subsurface)),
		title=title)
end

# ╔═╡ eb6f822d-9389-4cd9-adf2-75c29eb3b148
function plot_subsurface(ore; kwargs...)
	return plot_subsurface(zeros(size(ore)...), ore; kwargs...)
end

# ╔═╡ 586a1b2e-4c6a-4b5c-81ee-9f29ef3cd7a1
md"""
## Circles
"""

# ╔═╡ 0cf7c874-f583-4f89-8818-df478f25e142
@bind xy Slider([10, 20, 50, 80], default=50, show_value=true)

# ╔═╡ d6d62f08-ec32-45a6-943d-31057c033beb
σ_blur = xy / 25

# ╔═╡ d6106b64-2f8e-4f1d-914c-b8b0c94f41ac
dim_scale = 1/(xy/50)^2

# ╔═╡ a966f760-8c66-45dd-833a-b84e40644c1a
abstract type MainbodyGen end

# ╔═╡ 69fdb854-a8be-4e11-ab8e-be4ecf2824d5
abstract type ShapeNode <: MainbodyGen end

# ╔═╡ a28ad466-04fb-403e-b3e4-4917f25ab1c7
"""
Distribution for determining the center of the mainbody.
"""
function center_distribution(grid_dims; bounds=[grid_dims[1]/4, grid_dims[1]/2])
    xdistr = Distributions.Uniform(bounds[1], bounds[2])
    ydistr = Distributions.Uniform(bounds[1], bounds[2])
    return Product([xdistr, ydistr])
end

# ╔═╡ a6a7d16f-e34f-4f6a-a175-355f8c58d675
"""
Mainbody shape parameterized as a circle with `center` and `radius`.
"""
@with_kw struct CircleNode <: ShapeNode
    grid_dims::Tuple{Real, Real, Real}
    center::Union{Distribution, Vector} = center_distribution(grid_dims)
    radius::Union{Distribution, Real} = Distributions.Uniform(grid_dims[1]/12, grid_dims[1]/8)
end

# ╔═╡ d4873c89-7169-460c-be01-8c58e9825ea7
CircleNode <: ShapeNode <: MainbodyGen

# ╔═╡ 52a514d4-f96d-4cbb-85f2-556e58042c92
"""
Draw circle using Luxor. Relative to the translated `center`, hence `Point(0,0)`.
"""
function draw(shape::CircleNode)
    return circle(Luxor.Point(0,0), shape.radius, :fill)
end

# ╔═╡ 3244a070-5456-4097-a1aa-5b26a94460e9
Base.rand(shape::CircleNode; kwargs...) = rand(Random.GLOBAL_RNG, shape; kwargs...)

# ╔═╡ ee57dc48-a915-4277-88bb-74a5e92117bb
function clamp2dims!(x, dims)
	x[1] = clamp(x[1], 1, dims[1])
	x[2] = clamp(x[2], 1, dims[2])
	return x
end

# ╔═╡ 9498ae2d-dadc-452e-ab5e-cfd24957e7ad
md"""
## Perturb sample
"""

# ╔═╡ b00b0670-cf19-4776-8267-db965bddae00
@bind is_perturbed CheckBox()

# ╔═╡ b31d4955-5593-4a2c-b9f0-75cd974a8811
md"""
## Ellipse
"""

# ╔═╡ cc22c534-18c5-4d82-95b7-2ae63f9a40af
@with_kw struct EllipseNode <: ShapeNode
	grid_dims::Tuple{Real, Real, Real}
	center::Union{Distribution, Vector} = center_distribution(grid_dims)
	width::Union{Distribution, Real} = Distributions.Uniform(grid_dims[1]/6, grid_dims[1]/4)
	height::Union{Distribution, Real} = Distributions.Uniform(grid_dims[2]/6, grid_dims[2]/4)
	angle::Union{Distribution, Real} = Distributions.Uniform(0, 2π)
end

# ╔═╡ 0e914b68-2b92-4969-9003-5eaf3b03c000
function draw(shape::EllipseNode)
	return ellipse(Luxor.Point(0,0), shape.width, shape.height, :fill)
end

# ╔═╡ 19475725-9bac-415b-b428-a454f2d7d848
Base.rand(shape::EllipseNode; kwargs...) = rand(Random.GLOBAL_RNG, shape; kwargs...)

# ╔═╡ a285dbe2-03cf-4806-9373-42c4cbed1b5c
@bind is_perturbed_e CheckBox()

# ╔═╡ cfa3b14c-1726-4cca-b91e-cb6aa8a3da0d
md"""
## Bezier curve blobs
"""

# ╔═╡ dad99c6f-85db-47f5-a48e-324070dcd602
@with_kw struct BlobNode <: ShapeNode
	grid_dims::Tuple{Real, Real, Real}
	center::Union{Distribution, Vector} = center_distribution(grid_dims)
	N::Union{Distribution, AbstractArray, Real} = Distributions.Uniform(5, 8)
	factor::Union{Distribution, AbstractArray, Real} = Distributions.Uniform(5, 8)
	points::Union{Nothing, Vector{Luxor.Point}} = nothing
	angle::Union{Distribution, Real} = Distributions.Uniform(0, 2π)
end

# ╔═╡ 6609b356-2e9f-4e2a-a0a7-89ce7e7dd73b
function draw(shape::BlobNode)
	return drawbezierpath(makebezierpath(shape.points), :fill)
end

# ╔═╡ d3129f61-fe3f-42b4-879a-1e40586379ac
Base.rand(shape::BlobNode; kwargs...) = rand(Random.GLOBAL_RNG, shape; kwargs...)

# ╔═╡ 3e403f35-c513-4ecf-8a72-aa9e36742257
md"""
### Blob perturbation
"""

# ╔═╡ a9491ab5-3aaa-4282-94a1-c198641184d7
@bind is_perturbed_b CheckBox()

# ╔═╡ aa7801a7-d7e3-4a8f-bee6-cee34ff39a4f
md"""
## Multi-shape node
"""

# ╔═╡ bacb52cb-f10e-4f0f-a37a-5770f3524f59
struct MultiShapeNode <: MainbodyGen
	shapes::Vector{ShapeNode}
end

# ╔═╡ e172e1a2-3877-4c26-b901-065757683f84
node1 = EllipseNode(grid_dims=grid_dims)

# ╔═╡ 381f8808-70d7-4f62-bfea-d21715a80835
node2 = BlobNode(grid_dims=grid_dims, center=[35,35])

# ╔═╡ a4c6bd17-5baf-4d4b-bba9-82748364a90f
mainbody_m = MultiShapeNode([node1, node2])

# ╔═╡ da5c780e-f891-4b11-882e-f44722c5d211
function Base.rand(rng::Random.AbstractRNG, shape::MultiShapeNode; σ=2)
	shapes = []
	params = []
	for s in shape.shapes
		sampled_shape, sampled_params = rand(rng, s; σ=σ)
		push!(shapes, sampled_shape)
		push!(params, sampled_params)
	end
	return sum(shapes), params
end

# ╔═╡ bd9808a9-c072-4fcf-8167-aa7f6aed9934
Base.rand(shape::MultiShapeNode; kwargs...) = rand(Random.GLOBAL_RNG, shape; kwargs...)

# ╔═╡ 4c67eaa7-7b24-4af4-b332-7dbbda70c918
function perturb_sample(mainbody::MultiShapeNode, mainbody_params, noise)
	p_mainbody = []
	p_mainbody_params = []
	for (i,s) in enumerate(mainbody.shapes)
		p_mainbody_s, p_mainbody_params_s = perturb_sample(s, mainbody_params[i], noise)
		push!(p_mainbody, p_mainbody_s)
		push!(p_mainbody_params, p_mainbody_params_s)
	end
	
    return sum(p_mainbody), p_mainbody_params
end

# ╔═╡ 67acb8a9-a765-4780-980e-0965b9fd2c32
@bind is_perturbed_m CheckBox()

# ╔═╡ 10a098f2-de4a-4f80-a631-eba87c29eb7d
md"""
## Rectangle
"""

# ╔═╡ 0eaeb5fd-2664-49a2-a0f6-4d8d7b2b06a0
@with_kw struct RectangleNode <: ShapeNode
	grid_dims::Tuple{Real, Real, Real}
	center::Union{Distribution, Vector} = center_distribution(grid_dims)
	width::Union{Distribution, Real} = Distributions.Uniform(grid_dims[1]/6, grid_dims[1]/4)
	height::Union{Distribution, Real} = Distributions.Uniform(grid_dims[2]/6, grid_dims[2]/4)
	angle::Union{Distribution, Real} = Distributions.Uniform(0, 2π)
end

# ╔═╡ 32dab518-7bdd-49dd-a9a3-4dbf1b18ecb4
function draw(shape::RectangleNode)
	cornerpoint = Luxor.Point(-shape.width/2, -shape.height/2)
	return rect(cornerpoint, shape.width, shape.height, :fill)
end

# ╔═╡ 7c43a800-370f-42a2-8a10-c42f1ecbed40
Base.rand(shape::RectangleNode; kwargs...) = rand(Random.GLOBAL_RNG, shape; kwargs...)

# ╔═╡ 0c084275-f559-4ca9-9465-c5b234e765f2
@bind is_perturbed_r CheckBox()

# ╔═╡ e2f29813-29dc-475e-a514-c04becaa8a4b
md"""
## Generate shape matrix
"""

# ╔═╡ 78b10b13-e42b-4f18-990d-4158fadcc505
Base.convert(::Type{Float64}, m::ColorTypes.ARGB32) = convert(Float64, Gray(m))

# ╔═╡ b0c8e231-0b97-4538-9bca-3024d2dbd5de
blur(img, σ) = imfilter(img, Kernel.gaussian(σ))

# ╔═╡ 0bd7d2ef-5559-485d-9e82-33ca1eb53e03
function generate_shape_matrix(shape::ShapeNode; σ=2, grayscale=0.4)
    grid_dims = shape.grid_dims
    mat = @imagematrix begin
        gsave()
        background("black")
        Luxor.origin(Luxor.Point(0,0))
        pos = Luxor.Point(shape.center...)
        Luxor.translate(pos)
        if hasproperty(shape, :angle)
            Luxor.rotate(shape.angle)
        end
        sethue(grayscale, grayscale, grayscale)
        draw(shape)
        grestore()
    end grid_dims[1] grid_dims[2]

    return convert(Matrix{Float64}, Gray.(blur(mat, σ)))
end

# ╔═╡ 9a94006c-c0d4-4554-a475-649a9226716a
function Base.rand(rng::Random.AbstractRNG, shape::CircleNode; σ=2)
    grid_dims = shape.grid_dims
    center = isa(shape.center, Distribution) ? rand(rng, shape.center) : shape.center
    radius = isa(shape.radius, Distribution) ? rand(rng, shape.radius) : shape.radius
    shape = CircleNode(grid_dims=grid_dims, center=center, radius=radius)
    params = [center, radius, σ]
    return (generate_shape_matrix(shape; σ=σ), params)
end

# ╔═╡ 2dad36cb-6671-40f9-8e23-47986f3ab82f
function Base.rand(rng::Random.AbstractRNG, shape::EllipseNode; σ=2)
	grid_dims = shape.grid_dims
	center = isa(shape.center, Distribution) ? rand(shape.center) : shape.center
	width = isa(shape.width, Distribution) ? rand(shape.width) : shape.width
	height = isa(shape.height, Distribution) ? rand(shape.height) : shape.height
	angle = isa(shape.angle, Distribution) ? rand(shape.angle) : shape.angle
	shape = EllipseNode(grid_dims=grid_dims, center=center, width=width, height=height, angle=angle)
	params = [center, width, height, angle, σ]
	return (generate_shape_matrix(shape; σ=σ), params)
end

# ╔═╡ 0efc2ce4-6b6a-4328-b6d1-02bb03e1afe9
function Base.rand(rng::Random.AbstractRNG, shape::BlobNode; σ=2)
	grid_dims = shape.grid_dims
	center = isa(shape.center, Distribution) ? rand(shape.center) : shape.center
	N = isa(shape.N, Distribution) || isa(shape.N, AbstractArray) ? rand(shape.N) : shape.N
	factor = isa(shape.factor, Distribution) || isa(shape.factor, AbstractArray) ? rand(shape.factor) : shape.factor
	if isnothing(shape.points)
		pts = polysortbyangle(randompointarray(
			Luxor.Point(-grid_dims[1]/factor,-grid_dims[2]/factor),
			Luxor.Point(grid_dims[1]/factor, grid_dims[2]/factor),
			N))
	else
		pts = shape.points
	end
	angle = isa(shape.angle, Distribution) ? rand(shape.angle) : shape.angle
	shape = BlobNode(grid_dims=grid_dims, center=center, N=N, factor=factor, points=pts, angle=angle)
	params = [center, N, factor, pts, angle, σ]
	return (generate_shape_matrix(shape; σ=σ), params)
end

# ╔═╡ d9d694ea-fa58-48c6-9a12-535b8ae9fcc6
function Base.rand(rng::Random.AbstractRNG, shape::RectangleNode; σ=2)
	grid_dims = shape.grid_dims
	center = isa(shape.center, Distribution) ? rand(shape.center) : shape.center
	width = isa(shape.width, Distribution) ? rand(shape.width) : shape.width
	height = isa(shape.height, Distribution) ? rand(shape.height) : shape.height
	angle = isa(shape.angle, Distribution) ? rand(shape.angle) : shape.angle
	shape = RectangleNode(grid_dims=grid_dims, center=center, width=width, height=height, angle=angle)
	params = [center, width, height, angle, σ]
	return (generate_shape_matrix(shape; σ=σ), params)
end

# ╔═╡ da0416d0-1e6f-4135-87a7-e63c049417b2
begin
	Random.seed!(0)
	test_mainbody = CircleNode(grid_dims=(xy,xy,1))
	test_mainbody_map, test_mainbody_params = rand(test_mainbody; σ=σ_blur)
	plot_subsurface(test_mainbody_map, clims=true)
end

# ╔═╡ bfbf7967-3e2b-439e-9c96-1378e0eb489a
dim_scale * sum(test_mainbody_map .>= 0.3)

# ╔═╡ 8731eb7d-ee2b-48ee-ae8b-34c314c4897c
"""
Perturb shape by adding `noise` to its parameters.
"""
function perturb_sample(mainbody::CircleNode, mainbody_params, noise)
    grid_dims = mainbody.grid_dims
    center, radius, σ = mainbody_params
	noise_scale = grid_dims[1] / 50
    𝒟_noise = Distributions.Uniform(-noise_scale*noise, noise_scale*noise)

    p_center = center .+ rand(𝒟_noise, 2)
    clamp2dims!(p_center, grid_dims)
    p_radius = clamp(radius + rand(𝒟_noise), 0.1, Inf)

    p_shape = CircleNode(grid_dims=grid_dims, center=p_center, radius=p_radius)
    p_mainbody = generate_shape_matrix(p_shape; σ=σ)
    p_mainbody_params = [p_center, p_radius, σ]

    return p_mainbody, p_mainbody_params
end

# ╔═╡ 928fd048-6b0e-4342-9236-9c06a13b5900
begin
	Random.seed!(0)
	mainbody = CircleNode(grid_dims=(10,10,1), radius=0.5)
	mainbody_map, mainbody_params = rand(mainbody)
	plot_subsurface(mainbody_map)
end

# ╔═╡ 363d1af9-6e05-4d01-8052-dfed61816e6c
mainbody_map

# ╔═╡ d8610e13-b51b-4cdb-93fa-48903dc05db2
mainbody_params

# ╔═╡ 228c740f-b39e-44ca-8753-80fb056402b9
plot_subsurface(mean([rand(mainbody)[1] for _ in 1:500]), clims=false)

# ╔═╡ 2bbb8008-9e12-49ad-9d6c-4b0900646421
begin
	Random.seed!(0)
	mainbody_e = EllipseNode(grid_dims=grid_dims)
	mainbody_map_e, mainbody_params_e = rand(mainbody_e)
	plot_subsurface(mainbody_map_e)
end

# ╔═╡ bc1b94b9-df89-4f45-b1a0-5953f7e26636
function perturb_sample(mainbody::EllipseNode, mainbody_params, noise)
    grid_dims = mainbody.grid_dims
    center, width, height, angle, σ = mainbody_params
    𝒟_noise = Distributions.Uniform(-noise, noise)

    p_center = center .+ rand(𝒟_noise, 2)
    clamp2dims!(p_center, grid_dims)
    p_width = width + rand(𝒟_noise)
	p_height = height + rand(𝒟_noise)
	p_angle = angle + deg2rad(rand(𝒟_noise))

    p_shape = EllipseNode(grid_dims=grid_dims, center=p_center, width=p_width, height=p_height, angle=p_angle)
    p_mainbody = generate_shape_matrix(p_shape; σ=σ)
    p_mainbody_params = [p_center, p_width, p_height, p_angle, σ]

    return p_mainbody, p_mainbody_params
end

# ╔═╡ deb4191c-eb7c-411e-94ba-e0eb4a47b1e7
begin
	Random.seed!(0)
	mainbody_b = BlobNode(grid_dims=grid_dims)
	mainbody_map_b, mainbody_params_b = rand(mainbody_b; σ=2)
	plot_subsurface(mainbody_map_b)
end

# ╔═╡ 653dac11-a981-42a9-b973-f6e57901d9f9
mainbody_params_b

# ╔═╡ 03bd8438-3872-4011-af30-6795570a526f
points = mainbody_params_b[4]

# ╔═╡ 425e4fa0-dbfe-4d4b-b763-0d673445a0fa
[Luxor.Point(p.x + rand(), p.y + rand()) for p in points]

# ╔═╡ dff63bf1-90aa-41ab-b5ae-007051c00929
plot_subsurface(mean([rand(mainbody_b)[1] for _ in 1:500]), clims=false, title="blob")

# ╔═╡ ba4561da-021d-42e3-87ba-89604c135085
plot_subsurface(mean([rand(mainbody_e)[1] for _ in 1:500]), clims=false, title="ellipse")

# ╔═╡ 7ba5f080-a1f6-46b8-ab6e-bf74673bb201
plot_subsurface(mean([rand(mainbody)[1] for _ in 1:500]), clims=false, title="circle")

# ╔═╡ 1bb8e3a8-45ef-488c-879b-d979948f7ea8
function perturb_sample(mainbody::BlobNode, mainbody_params, noise)
    grid_dims = mainbody.grid_dims
    center, N, factor, points, angle, σ = mainbody_params
    𝒟_noise = Distributions.Uniform(-noise, noise)

    p_center = center .+ rand(𝒟_noise, 2)
    clamp2dims!(p_center, grid_dims)
    p_N = N + rand(𝒟_noise)
	p_factor = factor + rand(𝒟_noise)
	p_points = [Luxor.Point(p.x + rand(𝒟_noise), p.y + rand(𝒟_noise)) for p in points]
	p_angle = angle + deg2rad(rand(𝒟_noise))

    p_shape = BlobNode(grid_dims=grid_dims, center=p_center, N=p_N, factor=p_factor, points=p_points, angle=p_angle)
    p_mainbody = generate_shape_matrix(p_shape; σ=σ)
    p_mainbody_params = [p_center, p_N, p_factor, p_points, p_angle, σ]

    return p_mainbody, p_mainbody_params
end

# ╔═╡ e01a2b94-5fe6-45b5-85f3-4ebcfc4db494
begin
	Random.seed!(0)
	mainbody_map_m, mainbody_params_m = rand(mainbody_m; σ=0)
	plot_subsurface(mainbody_map_m)
end

# ╔═╡ 83e5735d-441c-4817-8ab4-7dc8b6f7bbe9
mainbody_params_m

# ╔═╡ 2cbc00a2-dd0e-49b7-8a36-9a7485442dda
begin
	Random.seed!(0)
	mainbody_r = RectangleNode(grid_dims=grid_dims, center=[20,10])
	mainbody_map_r, mainbody_params_r = rand(mainbody_r; σ=2)
	plot_subsurface(mainbody_map_r)
end

# ╔═╡ 16f5e689-aeb0-42b9-adf1-3a9bda37b4d2
mainbody_params_r

# ╔═╡ cb1933d2-7d1d-42d9-879e-39a2a34b0b44
function perturb_sample(mainbody::RectangleNode, mainbody_params, noise)
    grid_dims = mainbody.grid_dims
    center, width, height, angle, σ = mainbody_params
    𝒟_noise = Distributions.Uniform(-noise, noise)

    p_center = center .+ rand(𝒟_noise, 2)
    clamp2dims!(p_center, grid_dims)
    p_width = width + rand(𝒟_noise)
	p_height = height + rand(𝒟_noise)
	p_angle = angle + deg2rad(rand(𝒟_noise))

    p_shape = RectangleNode(grid_dims=grid_dims, center=p_center, width=p_width, height=p_height, angle=p_angle)
    p_mainbody = generate_shape_matrix(p_shape; σ=σ)
    p_mainbody_params = [p_center, p_width, p_height, p_angle, σ]

    return p_mainbody, p_mainbody_params
end

# ╔═╡ 5279b917-0754-47c2-ada8-08acd5823d5f
p_mainbody_map, p_mainbody_params = perturb_sample(mainbody, mainbody_params, 1.0);

# ╔═╡ 68c6a183-7b93-409f-a945-bad5cda23a86
p_mainbody_params

# ╔═╡ e1d03ae7-bf10-41fc-9de6-c29eb931afef
p_mainbody_map

# ╔═╡ 4cbb2ae7-82b8-43cc-abbf-befcb5df8c3a
if is_perturbed
	plot_subsurface(p_mainbody_map, clims=false)
	title!("perturbed")
else
	plot_subsurface(mainbody_map, clims=false)
	title!("original")
end

# ╔═╡ 754839af-5abd-4112-bea4-29ecf2d0ca6f
plot(
	begin
		plot_subsurface(p_mainbody_map, clims=false)
		title!("perturbed")
	end,
	begin
		plot_subsurface(mainbody_map, clims=false)
		title!("original")
	end,
	size=(600,250))

# ╔═╡ 9e7d661b-a8d0-4c48-962a-ee10bd1449c0
p_mainbody_map_e, p_mainbody_params_e = perturb_sample(mainbody_e, mainbody_params_e, 1.0);

# ╔═╡ fc91e2e9-0054-4564-a2e3-a05abdaca80e
if is_perturbed_e
	plot_subsurface(p_mainbody_map_e, clims=false)
	title!("perturbed")
else
	plot_subsurface(mainbody_map_e, clims=false)
	title!("original")
end

# ╔═╡ 0f94a2d7-59c6-4f8b-aae5-06f7c8af8383
p_mainbody_map_b, p_mainbody_params_b = perturb_sample(mainbody_b, mainbody_params_b, 1.0);

# ╔═╡ a8ad678c-716f-451d-9a22-53b00ca7fc9e
if is_perturbed_b
	plot_subsurface(p_mainbody_map_b, clims=false)
	title!("perturbed")
else
	plot_subsurface(mainbody_map_b, clims=false)
	title!("original")
end

# ╔═╡ 63f030e5-d1e1-41fb-9210-1f0483962fac
p_mainbody_map_m, p_mainbody_params_m = perturb_sample(mainbody_m, mainbody_params_m, 1.0);

# ╔═╡ 934a0548-e5f5-4b0d-ab35-a04641395168
if is_perturbed_m
	plot_subsurface(p_mainbody_map_m, clims=false)
	title!("perturbed")
else
	plot_subsurface(mainbody_map_m, clims=false)
	title!("original")
end

# ╔═╡ cbe219d6-4471-422e-ab3d-93d6c76d5ea4
p_mainbody_map_r, p_mainbody_params_r = perturb_sample(mainbody_r, mainbody_params_r, 1.0);

# ╔═╡ e4c0b6aa-e34c-480f-88f3-c7e055740d1f
if is_perturbed_r
	plot_subsurface(p_mainbody_map_r, clims=false)
	title!("perturbed")
else
	plot_subsurface(mainbody_map_r, clims=false)
	title!("original")
end

# ╔═╡ 23f8fbe9-a949-4c87-95fd-da07a2433470
function normalize01(mat)
	min = minimum(mat)
	max = maximum(mat)
	return (mat .- min) / (max - min)
end

# ╔═╡ a0b1b313-0ab0-4f35-98d1-4224152ee55f
md"""
---
"""

# ╔═╡ 4421d423-9c26-44ee-9866-d60a81ac95a6
# begin
# 	External = ingredients("..//..//MineralExploration//src//shapes.jl")
# 	import .External: TESTING
# end	

# ╔═╡ a5baa4a4-f0fc-4acf-a717-16dae5e92ad5
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
ImageFiltering = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
POMDPs = "a93abf59-7444-517b-a68a-c42f96afdd7d"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"

[compat]
Distributions = "~0.25.41"
ImageFiltering = "~0.7.1"
Luxor = "~3.0.0"
POMDPs = "~0.9.3"
Parameters = "~0.12.3"
Plots = "~1.25.7"
PlutoUI = "~0.7.32"
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

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

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

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "15bd05c1c0d5dbb32a9a3d7e0ad2d50dd6167189"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.1"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

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

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

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
# ╟─fdef4b60-7fbc-11ec-3adc-a7773dfccfa2
# ╠═9408fe18-99c9-48ff-9f6a-63caa2c99d96
# ╠═c2cfe5bc-ac18-4310-8a30-6b89c0a1e9d8
# ╠═539de988-9349-47b1-8129-d54d80663154
# ╠═e8b732da-65b3-4977-879d-87f256ec48c1
# ╠═0e1a8833-0ebd-4893-a463-7dc092a60540
# ╠═fc4e7fc0-9b2c-4e3c-80c9-7120ef1019c5
# ╟─8a71c40c-2319-479e-8fb9-710a57d540b3
# ╠═feb16f09-aba5-444b-a814-1c6a816a6c2f
# ╠═37fc09d9-cd50-4197-879a-045b6af2e707
# ╠═052eaac9-bac5-439a-8a56-058d81d7bb97
# ╠═6aba071e-5b44-43ae-b94c-1ad5bcd8baba
# ╠═e5bc2147-fcdf-449e-9693-3d902e27bbc1
# ╠═14fdb72c-28b4-4aab-b55b-32ed943b8ed8
# ╠═eb6f822d-9389-4cd9-adf2-75c29eb3b148
# ╟─586a1b2e-4c6a-4b5c-81ee-9f29ef3cd7a1
# ╠═0cf7c874-f583-4f89-8818-df478f25e142
# ╠═da0416d0-1e6f-4135-87a7-e63c049417b2
# ╠═d6d62f08-ec32-45a6-943d-31057c033beb
# ╠═d6106b64-2f8e-4f1d-914c-b8b0c94f41ac
# ╠═bfbf7967-3e2b-439e-9c96-1378e0eb489a
# ╠═a966f760-8c66-45dd-833a-b84e40644c1a
# ╠═69fdb854-a8be-4e11-ab8e-be4ecf2824d5
# ╠═d4873c89-7169-460c-be01-8c58e9825ea7
# ╠═a28ad466-04fb-403e-b3e4-4917f25ab1c7
# ╠═a6a7d16f-e34f-4f6a-a175-355f8c58d675
# ╠═52a514d4-f96d-4cbb-85f2-556e58042c92
# ╠═9a94006c-c0d4-4554-a475-649a9226716a
# ╠═3244a070-5456-4097-a1aa-5b26a94460e9
# ╠═ee57dc48-a915-4277-88bb-74a5e92117bb
# ╠═8731eb7d-ee2b-48ee-ae8b-34c314c4897c
# ╠═928fd048-6b0e-4342-9236-9c06a13b5900
# ╠═228c740f-b39e-44ca-8753-80fb056402b9
# ╠═363d1af9-6e05-4d01-8052-dfed61816e6c
# ╠═d8610e13-b51b-4cdb-93fa-48903dc05db2
# ╟─9498ae2d-dadc-452e-ab5e-cfd24957e7ad
# ╠═2baf8598-7964-4e77-9a6f-94427cf347bc
# ╠═5279b917-0754-47c2-ada8-08acd5823d5f
# ╠═68c6a183-7b93-409f-a945-bad5cda23a86
# ╠═e1d03ae7-bf10-41fc-9de6-c29eb931afef
# ╠═b00b0670-cf19-4776-8267-db965bddae00
# ╠═4cbb2ae7-82b8-43cc-abbf-befcb5df8c3a
# ╠═754839af-5abd-4112-bea4-29ecf2d0ca6f
# ╟─b31d4955-5593-4a2c-b9f0-75cd974a8811
# ╠═cc22c534-18c5-4d82-95b7-2ae63f9a40af
# ╠═0e914b68-2b92-4969-9003-5eaf3b03c000
# ╠═2dad36cb-6671-40f9-8e23-47986f3ab82f
# ╠═19475725-9bac-415b-b428-a454f2d7d848
# ╠═2bbb8008-9e12-49ad-9d6c-4b0900646421
# ╠═bc1b94b9-df89-4f45-b1a0-5953f7e26636
# ╠═9e7d661b-a8d0-4c48-962a-ee10bd1449c0
# ╠═a285dbe2-03cf-4806-9373-42c4cbed1b5c
# ╠═fc91e2e9-0054-4564-a2e3-a05abdaca80e
# ╟─cfa3b14c-1726-4cca-b91e-cb6aa8a3da0d
# ╠═dad99c6f-85db-47f5-a48e-324070dcd602
# ╠═6609b356-2e9f-4e2a-a0a7-89ce7e7dd73b
# ╠═0efc2ce4-6b6a-4328-b6d1-02bb03e1afe9
# ╠═d3129f61-fe3f-42b4-879a-1e40586379ac
# ╠═deb4191c-eb7c-411e-94ba-e0eb4a47b1e7
# ╠═653dac11-a981-42a9-b973-f6e57901d9f9
# ╠═03bd8438-3872-4011-af30-6795570a526f
# ╠═425e4fa0-dbfe-4d4b-b763-0d673445a0fa
# ╠═dff63bf1-90aa-41ab-b5ae-007051c00929
# ╠═ba4561da-021d-42e3-87ba-89604c135085
# ╠═7ba5f080-a1f6-46b8-ab6e-bf74673bb201
# ╟─3e403f35-c513-4ecf-8a72-aa9e36742257
# ╠═1bb8e3a8-45ef-488c-879b-d979948f7ea8
# ╠═0f94a2d7-59c6-4f8b-aae5-06f7c8af8383
# ╠═a9491ab5-3aaa-4282-94a1-c198641184d7
# ╠═a8ad678c-716f-451d-9a22-53b00ca7fc9e
# ╟─aa7801a7-d7e3-4a8f-bee6-cee34ff39a4f
# ╠═bacb52cb-f10e-4f0f-a37a-5770f3524f59
# ╠═e172e1a2-3877-4c26-b901-065757683f84
# ╠═381f8808-70d7-4f62-bfea-d21715a80835
# ╠═a4c6bd17-5baf-4d4b-bba9-82748364a90f
# ╠═da5c780e-f891-4b11-882e-f44722c5d211
# ╠═bd9808a9-c072-4fcf-8167-aa7f6aed9934
# ╠═e01a2b94-5fe6-45b5-85f3-4ebcfc4db494
# ╠═83e5735d-441c-4817-8ab4-7dc8b6f7bbe9
# ╠═4c67eaa7-7b24-4af4-b332-7dbbda70c918
# ╠═63f030e5-d1e1-41fb-9210-1f0483962fac
# ╠═67acb8a9-a765-4780-980e-0965b9fd2c32
# ╠═934a0548-e5f5-4b0d-ab35-a04641395168
# ╟─10a098f2-de4a-4f80-a631-eba87c29eb7d
# ╠═0eaeb5fd-2664-49a2-a0f6-4d8d7b2b06a0
# ╠═32dab518-7bdd-49dd-a9a3-4dbf1b18ecb4
# ╠═d9d694ea-fa58-48c6-9a12-535b8ae9fcc6
# ╠═7c43a800-370f-42a2-8a10-c42f1ecbed40
# ╠═2cbc00a2-dd0e-49b7-8a36-9a7485442dda
# ╠═16f5e689-aeb0-42b9-adf1-3a9bda37b4d2
# ╠═cb1933d2-7d1d-42d9-879e-39a2a34b0b44
# ╠═cbe219d6-4471-422e-ab3d-93d6c76d5ea4
# ╠═0c084275-f559-4ca9-9465-c5b234e765f2
# ╠═e4c0b6aa-e34c-480f-88f3-c7e055740d1f
# ╟─e2f29813-29dc-475e-a514-c04becaa8a4b
# ╠═0bd7d2ef-5559-485d-9e82-33ca1eb53e03
# ╠═78b10b13-e42b-4f18-990d-4158fadcc505
# ╠═b0c8e231-0b97-4538-9bca-3024d2dbd5de
# ╠═23f8fbe9-a949-4c87-95fd-da07a2433470
# ╟─a0b1b313-0ab0-4f35-98d1-4224152ee55f
# ╠═4421d423-9c26-44ee-9866-d60a81ac95a6
# ╟─a5baa4a4-f0fc-4acf-a717-16dae5e92ad5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
