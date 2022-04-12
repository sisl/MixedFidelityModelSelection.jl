### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 8405e36f-b6cf-481c-bb51-b67b67372c92
using StatsPlots

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

# ╔═╡ 41df3571-92bb-4520-b021-73dfd695d07d
using DataFrames

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

# ╔═╡ ec1e38f4-6e90-4041-8a96-b1d829de193c
results = BSON.load("..\\scripts\\MEParallel.jl\\results\\results_fixed_bores.bson")[:results]

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

# ╔═╡ 3085e225-59b9-4f0c-9586-5d82f307bf34
md"""
## Ore maps
"""

# ╔═╡ b271a899-2831-456d-b5ad-5868a5890bab
ore_map = mean(results50x50blob[:ore_map]);

# ╔═╡ f0f3efb8-e52d-4552-9675-762c4f3fedc4
begin
	mass_map = mean(results50x50blob[:mass_map])
	r_massive = mean(results50x50blob[:r_massive])
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

# ╔═╡ 4e4f6167-875a-41a7-b84a-53d65f5c00a7
md"""
## Relative volume errors
"""

# ╔═╡ 620df0d7-3034-4d90-b199-84f323903f22
dig2(x) = round(x, digits=2)

# ╔═╡ b111b064-56a0-42fe-8379-155dadd80f04
smry = x->(dig2(mean(x)),dig2(std(x)))

# ╔═╡ e00d2a2e-2258-4121-b96d-a80e4f30d7e1
begin

μb50, σb50 = smry(last.(results50x50blob[:rel_errors]))
μe50, σe50 = smry(last.(results50x50ellipse[:rel_errors]))
μc50, σc50 = smry(last.(results50x50circle[:rel_errors]))
μb30, σb30 = smry(last.(results30x30blob[:rel_errors]))
μe30, σe30 = smry(last.(results30x30ellipse[:rel_errors]))
μc30, σc30 = smry(last.(results30x30circle[:rel_errors]))

sb50 = "\$$μb50 \\pm $σb50\$"
se50 = "\$$μe50 \\pm $σe50\$"
sc50 = "\$$μc50 \\pm $σc50\$"
sb30 = "\$$μb30 \\pm $σb30\$"
se30 = "\$$μe30 \\pm $σe30\$"
sc30 = "\$$μc30 \\pm $σc30\$"

Markdown.parse("""

 –    | blob  | ellipse | circle
:---- | :---- | :------ | :-----
50×50 | $sb50 | $se50   | $sc50
30×30 | $sb30 | $se30   | $sc30
""")

end

# ╔═╡ 93545762-a44d-41a8-8e00-e1493d0748e4
plot_rel_error_aggregate(results50x50blob, results50x50ellipse, results50x50circle)

# ╔═╡ 378681c1-8754-4381-bf94-fa956dc58e70
plot_rel_error_aggregate(results30x30blob, results30x30ellipse, results30x30circle)

# ╔═╡ 45ad5118-df2e-44b9-9a24-20067944b06c
md"""
## Polar comparison plots
- accuracy
- runtime
- regret
"""

# ╔═╡ ce02b5d0-4f40-4dcc-897d-8f49082725af
md"""
## Confusion statistics
"""

# ╔═╡ 80afea6c-68b1-4f4d-96c6-bae16f3bc8ac
decisions = results50x50blob[:last_action]

# ╔═╡ 616e0bb0-3944-4f33-9216-afd40015382d
true_decisions = get_true_decisions(results50x50blob)

# ╔═╡ 7229656e-3525-4be4-a11f-dd277c3b46cd
cm = confusion_matrix(decisions, true_decisions)

# ╔═╡ 02c3f953-69b7-4419-96d2-18341d8f47f8
plot_confusion(cm)

# ╔═╡ f2be7470-82d6-4133-821a-68c887a4a2f0
begin
	dtd = res->(res[:last_action], get_true_decisions(res))
	
	db50, tdb50 = dtd(results50x50blob)
	de50, tde50 = dtd(results50x50ellipse)
	dc50, tdc50 = dtd(results50x50circle)
	db30, tdb30 = dtd(results30x30blob)
	de30, tde30 = dtd(results30x30ellipse)
	dc30, tdc30 = dtd(results30x30circle)
	
	cmb50 = confusion_matrix(db50, tdb50)
	cme50 = confusion_matrix(de50, tde50)
	cmc50 = confusion_matrix(dc50, tdc50)
	cmb30 = confusion_matrix(db30, tdb30)
	cme30 = confusion_matrix(de30, tde30)
	cmc30 = confusion_matrix(dc30, tdc30)
end;

# ╔═╡ 60097607-f34d-4a95-a843-efeac51297a3
confusion_table(cmb50, :blob, 50)

# ╔═╡ eccd50b2-4401-4cde-b2b6-021ca5263878
confusion_table(cme50, :ellipse, 50)

# ╔═╡ 0326ec74-2f6d-4d9b-93e1-fa5a0f3422a8
confusion_table(cmc50, :circle, 50)

# ╔═╡ 29e760d4-3e7d-49ca-9698-d9c12619f720
confusion_table(cmb30, :blob, 30)

# ╔═╡ 91c09551-3a27-4fcc-b5d7-811a652fc7db
confusion_table(cme30, :ellipse, 30)

# ╔═╡ 54cd37cf-9e53-43f5-bc3e-237cd7942c00
confusion_table(cmc30, :circle, 30)

# ╔═╡ e9854a55-31ad-4d7f-8331-e5575a6375c6
md"""
# 3D value plot
"""

# ╔═╡ ce739257-3473-426f-aef4-86a7e016a575
md"""
## Bayesian analysis using beta distributions
"""

# ╔═╡ 44df4b43-e782-47fa-8c64-f48162aa4bd8
results_regret = BSON.load("..\\scripts\\MEParallel.jl\\results\\results_10K_blobbiasfix.bson")[:results]

# ╔═╡ cf562730-8ac6-4b45-a311-a8208c3982fb
shapekeys = [:blob, :ellipse, :circle]

# ╔═╡ 83a79501-d3ba-4adf-b183-4d1b28b815bd
md"""
## Plotting accuracies
"""

# ╔═╡ a3ba5f18-7d8e-4a57-85e1-56933771bd77
plot_accuracy(results_regret, :blob)

# ╔═╡ 9e716e52-b402-442d-9ee4-36612d2d72d0
plot_accuracies(results_regret, shapekeys)

# ╔═╡ 53d2a31e-0fd5-49a9-ae81-5a55b22e328d
md"""
## Regret

$$\begin{gather}
R_\text{best} = \max\left\{0, \text{massive ore} - \text{extraction cost}\right\}\\
\operatorname{regret} = R_\text{best} - R
\end{gather}$$
"""

# ╔═╡ 157c59f3-8034-4575-9e96-3ef3025170be
begin
	plot_cumulative_regret(results_regret, (:ellipse, (10,10,1), 10); α_quantile=0.9)
	plot!(size=(350,200))
	xlims!(0, 150)
end

# ╔═╡ 74fdc433-91ff-4c08-9c46-2563be98f317
regret_fn_90th_quant = (res)->quantile(regret(res), 0.9)

# ╔═╡ b3296ce7-f606-4c43-b10a-3603c887be4e
begin
	regret_fn_mean = res->mean(regret(res))
	regret_title_mean = "expected regret"
	p_regret = plot_sweep_regret(results_regret, shapekeys, regret_fn_mean,
		regret_title_mean; colorbar_regret_fn=regret_fn_90th_quant)
end

# ╔═╡ eaebf4df-c915-4136-aa99-7101563d4aca
p_regret_color = plot_sweep_regret(results_regret, shapekeys, regret_fn_mean, regret_title_mean)

# ╔═╡ 94b106f6-409d-457f-a4a7-4ad68a9fd0ef
begin
	α_quantile = 0.8
	regret_fn_quant = res->quantile(regret(res), α_quantile)
	regret_title_quant = "$(round(Int, 100*α_quantile))th percentile of regret"
	p_regret80 = plot_sweep_regret(results_regret, shapekeys, regret_fn_quant, 
		regret_title_quant; colorbar_regret_fn=regret_fn_90th_quant)
end

# ╔═╡ 3bad29de-edc2-4248-b5fb-5c1aa360c457
begin
	α_quantile2 = 0.9
	regret_fn_quant90 = res->quantile(regret(res), α_quantile2)
	regret_title_quant90 = "$(round(Int, 100*α_quantile2))th percentile of regret"
	p_regret90 = plot_sweep_regret(results_regret, shapekeys, regret_fn_quant90, 
		regret_title_quant90; colorbar_regret_fn=regret_fn_90th_quant)
end

# ╔═╡ c1bccfa7-9f45-4ed2-8a4b-5a269329ed0a
quantile(regret(results_regret[:blob, (50,50,1), 10]), 0.2)

# ╔═╡ 9b3036a7-26ac-459f-9150-3cf537258744
plot_regret(results_regret, :blob)

# ╔═╡ c21fe459-c5c6-4324-ac17-005dddc94bb7
md"""
# Pareto curves (regret vs. runtime)
"""

# ╔═╡ e92ed5d4-955f-4138-88df-878448f1bd72
p_pareto, pareto_optimal = plot_pareto(results_regret; minutes=true, return_optimal=true);

# ╔═╡ df4bdef7-6a67-4933-b2a5-121feff6c0ee
p_pareto

# ╔═╡ 0a2c5707-6b04-48f9-ac2d-200392a62752
pareto_optimal

# ╔═╡ d76fd9c4-bc40-47f9-b6c5-7334e5d2a997
p_pareto_sec, pareto_optimal_sec = plot_pareto(results_regret; minutes=false, return_optimal=true);

# ╔═╡ cd15ed2f-a280-4bf3-9d7e-3bb9afa4720d
p_pareto_sec

# ╔═╡ 49cfde98-23b4-48d0-bed3-d739ad3b79f2
pareto_optimal_sec

# ╔═╡ 1d2062bc-63ee-4103-b476-e1b7aab6c453
md"""
# Ideas
"""

# ╔═╡ bc1948a9-b5f2-429d-8078-c1baca0e3482
p_regret_color

# ╔═╡ 98c42154-51c7-4274-bd46-3d53166f0575
p_regret90

# ╔═╡ 34bec962-d352-4bf1-a73d-993943def892
struct FidelityLevel
	level
	grid_dims
	pomcpow_iters
end

# ╔═╡ 7344ad53-3823-4a09-82b3-2600b78f66d8
fidelities = [[(10,10,1), (30,30,1), (50,50,1)], [10, 100, 1000]]

# ╔═╡ 0c48c759-9203-47ee-8442-74ab7b3d5746
shape = :ellipse

# ╔═╡ 1b92eaa2-63c8-4cb9-9819-05284eda89fd
begin
	F = []
	for (i,gd) in enumerate(fidelities[1])
		for (j,iter) in enumerate(fidelities[2])
			push!(F, FidelityLevel(i*j, gd, iter))
		end
	end
	F
end

# ╔═╡ c8c7b01a-a43c-49a9-b2fd-6507b97d4130
function create_dataframe(data, labels)
	series = vcat([fill(labels[i], length(data[i])) for i in 1:length(labels)]...)
	return DataFrame(series=series, data=vcat(data...))
end

# ╔═╡ 44dda21d-0fdd-438e-aefc-52d25397b3f2
begin
	D1 = regret(results_regret[(shape, (10,10,1), 10)])
	D2 = regret(results_regret[(shape, (10,10,1), 100)])
	D3 = regret(results_regret[(shape, (10,10,1), 1000)])
	df = create_dataframe([D1, D2, D3], [:D1, :D2, :D3])
	@df df groupedhist(:data, group=:series)	
end

# ╔═╡ 82bbf51b-347e-4b8d-b713-c87779eb7b3d
md"""
# Aggregate results
"""

# ╔═╡ 97cfc22a-e2fc-4ea4-b8a5-ffbeb3f4615a
function getall(results; shape=nothing, grid_dims=nothing, pomcpow_iters=nothing)
	res = Dict()
	for k in keys(results)
		candidate = nothing
		if !isnothing(shape) && k[1] == shape
			candidate = results[k]
		end
		if !isnothing(grid_dims) && k[2] != grid_dims
			candidate = nothing
		end
		if !isnothing(pomcpow_iters) && k[3] != pomcpow_iters
			candidate = nothing
		end
		if !isnothing(candidate)
			res[k] = candidate
		end
	end
	return res
end

# ╔═╡ 3108dcb0-c1cb-4597-8e23-8652fede05a7
begin
	res_c_10 = getall(results_regret; shape=:circle, pomcpow_iters=10)
	res_c_100 = getall(results_regret; shape=:circle, pomcpow_iters=100)
	res_c_1000 = getall(results_regret; shape=:circle, pomcpow_iters=1000)
	
	res_e_10 = getall(results_regret; shape=:ellipse, pomcpow_iters=10)
	res_e_100 = getall(results_regret; shape=:ellipse, pomcpow_iters=100)
	res_e_1000 = getall(results_regret; shape=:ellipse, pomcpow_iters=1000)

	res_b_10 = getall(results_regret; shape=:blob, pomcpow_iters=10)
	res_b_100 = getall(results_regret; shape=:blob, pomcpow_iters=100)
	res_b_1000 = getall(results_regret; shape=:blob, pomcpow_iters=1000)
end;

# ╔═╡ 1a5cd3f9-ad73-49e3-baf2-cf82e1a5726f
begin
    res_c_10x = getall(results_regret; shape=:circle, grid_dims=(10,10,1))
    res_c_30x = getall(results_regret; shape=:circle, grid_dims=(30,30,1))
    res_c_50x = getall(results_regret; shape=:circle, grid_dims=(50,50,1))
    
    res_e_10x = getall(results_regret; shape=:ellipse, grid_dims=(10,10,1))
    res_e_30x = getall(results_regret; shape=:ellipse, grid_dims=(30,30,1))
    res_e_50x = getall(results_regret; shape=:ellipse, grid_dims=(50,50,1))

    res_b_10x = getall(results_regret; shape=:blob, grid_dims=(10,10,1))
    res_b_30x = getall(results_regret; shape=:blob, grid_dims=(30,30,1))
    res_b_50x = getall(results_regret; shape=:blob, grid_dims=(50,50,1))
end;

# ╔═╡ c2110562-836f-4e4d-a999-b0fffe211da0
μσ(x) = (mean(x), std(x))

# ╔═╡ 8813ffa1-9a48-47c1-a096-5f88dd9b4a7c
[μσ(last.(res[:rel_errors])) for (k,res) in res_c_1000]

# ╔═╡ 0ae34129-c939-4d2a-8e52-519c806ff5a5
[(mean(regret(res)), std(regret(res))) for res in values(res_c_10)]

# ╔═╡ d3542175-9035-4dba-b2dd-257848e07510
function plot_rel_error_pomcpow(res; offset=0, c=:crimson, m=:circle)
	err10 = μσ(last.(merge(values(res[1])...)[:rel_errors]))
	err100 = μσ(last.(merge(values(res[2])...)[:rel_errors]))
	err1000 = μσ(last.(merge(values(res[3])...)[:rel_errors]))
	scatter!([10+offset], [err10[1]], label=false, c=c, marker=m, yerr=err10[2])
	scatter!([100+10*offset], [err100[1]], label=false, c=c, marker=m, yerr=err100[2])
	scatter!([1000+100*offset], [err1000[1]], label=false, c=c, marker=m, yerr=err1000[2])
	plot!(xaxis=:log)
end

# ╔═╡ 62e22eea-a560-4ba0-a69a-cea65301e669
function plot_rel_error_pomcpow2(res; offset=0, c=:crimson, m=:circle)
	err10 = μσ(last.(merge(values(res[1])...)[:rel_errors]))
	err100 = μσ(last.(merge(values(res[2])...)[:rel_errors]))
	err1000 = μσ(last.(merge(values(res[3])...)[:rel_errors]))
	plot!([10+offset, 100+10*offset, 1000+100*offset], [err10[1], err100[1], err1000[1]], label=false, c=c, marker=m, yerr=[err10[2], err100[2], err1000[2]])
	# scatter!([], [err100[1]], label=false, c=c, marker=m, yerr=err100[2])
	# scatter!([1000+100*offset], [err1000[1]], label=false, c=c, marker=m, yerr=err1000[2])
	plot!(xaxis=:log)
end

# ╔═╡ c66a8213-6d50-4bfa-b9d3-6e0e03b8b948
begin
	plot()
	plot_rel_error_pomcpow2([res_c_10, res_c_100, res_c_1000]; offset=-1, m=:circle)
	plot_rel_error_pomcpow2([res_e_10, res_e_100, res_e_1000]; offset=0, m=:diamond, c=:forestgreen)
	plot_rel_error_pomcpow2([res_b_10, res_b_100, res_b_1000]; offset=1, m=:star6, c=:blue)
    annotate!(15, 20, text("overestimated", :gray, :left, 10, "Computer Modern"))
    annotate!(15, -20, text("underestimated", :gray, :left, 10, "Computer Modern"))
	# draw_zero()
	plot!([xlims()...], [0, 0], lw=1, c=:gray, xlims=xlims(), label=false, style=:dash)
	xlabel!("planning iterations")
	ylabel!("relative error (at end of episode)")
	title!("relative error over planning iterations")
	ylims!(ylims()[1]*2, ylims()[2]*2)
end

# ╔═╡ 932679df-b6e9-416d-89bb-f9e130e6a1f7
function plot_rel_error_grid(res; offset=0, c=:crimson, m=:circle)
	err10 = μσ(last.(merge(values(res[1])...)[:rel_errors]))
	err30 = μσ(last.(merge(values(res[2])...)[:rel_errors]))
	err50 = μσ(last.(merge(values(res[3])...)[:rel_errors]))
	plot!([10+offset, 30+offset, 50+offset], [err10[1], err30[1], err50[1]], label=false, c=c, marker=m, yerr=[err10[2], err30[2], err50[2]])
	# scatter!([], [err100[1]], label=false, c=c, marker=m, yerr=err100[2])
	# scatter!([1000+100*offset], [err1000[1]], label=false, c=c, marker=m, yerr=err1000[2])
end

# ╔═╡ e784ecf7-d5b4-4b57-93db-db019f69671a
begin
    plot()
    plot_rel_error_grid([res_c_10x, res_c_30x, res_c_50x]; offset=-1, m=:circle)
    plot_rel_error_grid([res_e_10x, res_e_30x, res_e_50x]; offset=0, m=:diamond, c=:forestgreen)
    plot_rel_error_grid([res_b_10x, res_b_30x, res_b_50x]; offset=1, m=:star6, c=:blue)
    annotate!(15, 20, text("overestimated", :gray, :left, 10, "Computer Modern"))
    annotate!(15, -20, text("underestimated", :gray, :left, 10, "Computer Modern"))
    plot!([xlims()...], [0, 0], lw=1, c=:gray, xlims=xlims(), label=false, style=:dash)
    xlabel!("grid dimensions")
    ylabel!("relative error (at end of episode)")
    title!("relative error over grid dimensions")
    ylims!(ylims()[1]*2, ylims()[2]*2)
end

# ╔═╡ c95a26dd-2600-4a0e-acf2-9804cac84e8b
μσ(regret(merge(values(res_c_10x)...)))

# ╔═╡ fce48123-21bc-4312-a76b-86592d436dc4
function plot_regret_grid(res; offset=0, c=:crimson, m=:circle)
    regret10 = μσ(regret(merge(values(res[1])...)))
    regret30 = μσ(regret(merge(values(res[2])...)))
    regret50 = μσ(regret(merge(values(res[3])...)))
    plot!([10+offset, 30+offset, 50+offset], [regret10[1], regret30[1], regret50[1]], label=false, c=c, marker=m, yerr=[regret10[2], regret30[2], regret50[2]])
end

# ╔═╡ 34e518bb-6c7c-4dff-ab22-7356dfd5f6fb
begin
    plot()
    plot_regret_grid([res_c_10x, res_c_30x, res_c_50x]; offset=-1, m=:circle)
    plot_regret_grid([res_e_10x, res_e_30x, res_e_50x]; offset=0, m=:diamond, c=:forestgreen)
    plot_regret_grid([res_b_10x, res_b_30x, res_b_50x]; offset=1, m=:star6, c=:blue)
    xlabel!("grid dimensions")
    ylabel!("regret")
    title!("regret over grid dimensions")
    ylims!(ylims()[1]*2, ylims()[2]*2)
end

# ╔═╡ 61252668-bac4-4fdc-9cef-6f00fe7cb235
function plot_regret_pomcpow(res; offset=0, c=:crimson, m=:circle)
    regret10 = μσ(regret(merge(values(res[1])...)))
    regret100 = μσ(regret(merge(values(res[2])...)))
    regret1000 = μσ(regret(merge(values(res[3])...)))
    plot!([10+offset, 100+10*offset, 1000+100*offset], [regret10[1], regret100[1], regret1000[1]], label=false, c=c, marker=m, yerr=[regret10[2], regret100[2], regret1000[2]])
    plot!(xaxis=:log)
end

# ╔═╡ f09b7aa9-d48a-465c-a246-094f168283fd
begin
    plot()
    plot_regret_pomcpow([res_c_10, res_c_100, res_c_1000]; offset=-1, m=:circle)
    plot_regret_pomcpow([res_e_10, res_e_100, res_e_1000]; offset=0, m=:diamond, c=:forestgreen)
    plot_regret_pomcpow([res_b_10, res_b_100, res_b_1000]; offset=1, m=:star6, c=:blue)
    xlabel!("planning iterations")
    ylabel!("regret")
    title!("regret over planning iterations")
    ylims!(ylims()[1]*2, ylims()[2]*2)
end

# ╔═╡ a3c97203-6bf6-4903-ae26-eb2ed326ed08
MFMS = MixedFidelityModelSelection

# ╔═╡ 55ef731c-fab2-4371-9eea-f6df5cb16393
begin
	fig = MFMS.Figure(font="Computer Modern", resolution=(600,400))
	latexize = str -> MFMS.latexstring("\\mathrm{$(replace(str, " "=>"\\; "))}")

	offset_scale = 6

	ylab = "grid dimensions"
	shapes = latexize.(["10x10", "30x30", "50x50"])
	# fulldata = [res_c_10x, res_c_30x, res_c_50x]
 #    titleshape = "circle"
    
    # fulldata = [res_e_10x, res_e_30x, res_e_50x]
    # titleshape = "ellipse"
    
    fulldata = [res_b_10x, res_b_30x, res_b_50x]
    titleshape = "blob"
	
	# ylab = "planning iterations"
    # shapes = latexize.(["10", "100", "1000"])
	# fulldata = [res_c_10, res_c_100, res_c_1000]
	# titleshape = "circle"
	
	# fulldata = [res_e_10, res_e_100, res_e_1000]
	# titleshape = "ellipse"
	
	# fulldata = [res_b_10, res_b_100, res_b_1000]
	# titleshape = "blob"

	axes = MFMS.Axis(fig[1,1], title=latexize("regret distributions ($titleshape belief)"),
		 yticks=((1:length(shapes)) ./ offset_scale, reverse(shapes)),
		 xlabel=latexize("regret"),
		 ylabel=latexize(ylab),
		 xtickformat = xs -> [latexize(string(Int(x))) for x in xs])

	for i in 1:length(fulldata)
		data = convert(Vector{Real}, regret(merge(values(reverse(fulldata)[i])...)))
		
		MFMS.density!(data, offset=i/offset_scale, color=:x, colormap=:viridis, strokewidth=1, strokecolor=:black)

		μ = mean(data)
        σ = std(data)
        Xs = [μ, μ]
        Ys = [i/offset_scale, (i+0.25)/offset_scale]
        # Ys[2] = i == 3 ? Ys[2]*1.05 : Ys[2]
        MFMS.lines!(Xs, Ys, color=:red, linewidth=3)
        MFMS.text!(MFMS.latexstring("$(round(μ, digits=2)) ± $(round(σ, digits=2))"),
                         position=(Xs[2], Ys[2]+(0.05/offset_scale)),
                         align=(:left, :baseline),
                         textsize=14, color=:black)
	end
			
end

# ╔═╡ 75678d11-d516-48b5-af1b-2b31e4c5ce83
fig

# ╔═╡ 81eb626e-3e37-47d1-b64e-8cb9e5adc9fd
md"""
# Bias and variance
"""

# ╔═╡ 09533724-55ea-46a5-b016-3b558d452c07
f(res, extraction_cost=150) = res[:r_massive]

# ╔═╡ a924217b-127a-47e0-9f4e-af4e38f40ecb
f̂(res) = f(res) .+ last.(res[:rel_errors])

# ╔═╡ ec47b44d-5d53-4081-8cf7-a9cf0e923c88
# \$\\operatorname{bias}[\\hat{f}(x)] = \\mathbb{E}[\\hat{f}(x)] - f(x)\$
"""
\$\\operatorname{bias}[\\hat{f}(x)] = \\mathbb{E}[\\mathbb{E}[\\hat{f}(x)] - f(x)]\$
"""
bias(res) = mean(mean(f̂(res)) .- f(res))

# ╔═╡ 749ced45-e465-407f-b4c1-733e03802e6b
"""
\$\\operatorname{var}[\\hat{f}(x)] = \\mathbb{E}[(\\hat{f}(x) - \\mathbb{E}[\\hat{f}(x)])^2]\$
"""
variance(res) = mean((f̂(res) .- mean(f̂(res))).^2)

# ╔═╡ ca17880a-ceca-42a4-82e4-6ccc16590278
"""
\$\\operatorname{MSE} = \\mathbb{E}[(y - \\mathbb{E}[\\hat{f}(x)])^2]\$
"""
mse(res) = mean((f(res) - f̂(res)).^2)

# ╔═╡ f97c751a-5ebd-4233-85a6-5a7dc59dcd09
stddev(res) = sqrt(variance(res))

# ╔═╡ 189c3d48-7f25-43d6-b10d-b51493b8af7e
bores(res) = mean(res[:n_drills])

# ╔═╡ 65d715a5-2124-4388-b901-8419a0e786a4
returns(res) = mean(res[:discounted_return])

# ╔═╡ 33db8806-1e0c-4257-80f8-a4e4a5f2cff5
returns_var(res) = var(res[:discounted_return])

# ╔═╡ eb73d59a-631d-4bca-8f94-a47f6f6137ef
res = results_regret[(:blob, (50, 50, 1), 1000)]

# ╔═╡ f60da90b-daf6-49ff-991a-4e7233d3b56b
f(res)

# ╔═╡ fea21277-d833-4828-b980-70f550776a01
f̂(res)

# ╔═╡ 5232cdc7-5947-4dd7-a054-955f303983ca
bias(res)

# ╔═╡ 14a7c259-f243-4242-8913-34deba5c0659
stddev(res)

# ╔═╡ 04cac42a-a4f9-4c0f-b218-bfe199c4a955
function plot_fidelity_value(results, shapekey; value_func, kwargs...)
    plot()

    X = sort(unique(map(k->k[2][1], collect(keys(results))))) # (x,y,z) grab x
    Y = sort(unique(map(k->k[3], collect(keys(results))))) # POMCPOW iterations

    contourf!(X, Y, (x,y)->value_func(results[(shapekey, (x,x,1), y)]);
              c=:viridis, levels=16, size=(500,400), yaxis=:log, linewidth=0.25,
			  linecolor=:black, kwargs...)

    saved_lims = (xlims(), ylims())
    G = [(x, y) for x in X, y in Y]
    _X, _Y = first.(G), last.(G)
    scatter!(_X, _Y, ms=4, color=:white, label=false, marker=:circle)

    return plot!()
end

# ╔═╡ 0cd25377-6711-47b3-9c33-35d7df0f0640
begin
plot_biases(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, bias; title="bias", kwargs...)
plot_variances(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, variance; title="variance", kwargs...)
plot_stddev(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, stddev; title="standard deviation", kwargs...)
plot_mse(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, mse; title="mean squared error (MSE)", kwargs...)
plot_bores(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, bores; title="num. bores", kwargs...)
plot_returns(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, returns; title="discounted return", kwargs...)
plot_returns_var(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, returns_var; title="discounted return (variance)", kwargs...)

function plot_fidelities(results, shapekeys, value_fn; colorbar_fn=value_fn, title="bias", kwargs...)
    shared_cmap, svmin, svmax = get_colorbar_bounded_multi(results; value_fn=colorbar_fn)

    if value_fn == bias
        shared_cmap = get_colorbar_bias(svmin, svmax, vmid=0) # :broc # :curl
	elseif value_fn == variance || value_fn == stddev || value_fn == returns_var
        shared_cmap = get_colorbar_var(svmin, svmax)
    elseif value_fn == mse
        shared_cmap = get_colorbar_mse(svmin, svmax)
    elseif value_fn == bores
        shared_cmap = get_colorbar_bores(svmin, svmax)
    elseif value_fn == returns
        shared_cmap = get_colorbar_returns(svmin, svmax)
    end

    clims = (svmin, svmax)
    pbfunc = (s; kwargs...) -> begin
        plot_fidelity(results, shapekeys[s]; reduced=true, value_func=value_fn,
                  cmap=shared_cmap, clims=clims, kwargs...)
    end

    ptitle = plot(title="$title", grid=false, axis=false, tick=nothing, top_margin=-5Plots.mm)

    pplot = plot([pbfunc(s; show_ylabel=(s==1), show_xlabel=(s==2))
                 for s in 1:length(shapekeys)]...,
        layout=@layout[A B C])
    
    pcbar = contourf([0], [0], (x,y)->0,
              clims=clims, levels=15,
              c=shared_cmap, axis=false, tick=nothing, label=false)
    # pcbar = scatter([0], [0], alpha=0,
    #     zcolor=[clims...], clims=clims, c=shared_cmap,
    #     axis=false, tick=nothing, label=false)
    plot(ptitle, pplot, pcbar,
         layout=@layout([a{0.01h}; b c{0.1w}]),
         size=(710,250), bottom_margin=5Plots.mm)
end


function plot_fidelity(results, shapekey;
        reduced=false, cmap=nothing, clims=nothing,
        show_cbar=!reduced, show_xlabel=!reduced, show_ylabel=!reduced,
        value_func=res->bias(res))
    plot()
    title = "$shapekey"
    if reduced
        title!(title)
    else
        title!("? [$title]")
    end

    if show_xlabel
        xlabel!("grid dimension (n×n)")
    end
    if show_ylabel
        ylabel!("planning iterations")
    end

    X, Y = get_data_xy(results)

    if isnothing(cmap)
        # Based off individual values of the data (i.e., not shared)
        cmap = :broc
    end

    contourf!(X, Y, (x,y)->value_func(results[(shapekey, (x,x,1), y)]),
              c=cmap, cbar=show_cbar, clims=clims, levels=15,
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

function get_colorbar_bias(vmin, vmax; vmid=0)
    colors = [:black, :saddlebrown, :white, :cornflowerblue, :darkblue]
    return get_colorbar_normalized(vmin, vmax; colors=colors, vmid=vmid)
end

function get_colorbar_var(vmin, vmax)
    colors = :YlOrBr # [:forestgreen, :lightgreen, :gold, :lightcoral, :darkred]
    return get_colorbar_normalized(vmin, vmax; colors=colors)
end

function get_colorbar_mse(vmin, vmax)
    colors = :Reds
    return get_colorbar_normalized(vmin, vmax; colors=colors)
end

function get_colorbar_bores(vmin, vmax)
    colors = :viridis
    return get_colorbar_normalized(vmin, vmax; colors=colors)
end

function get_colorbar_returns(vmin, vmax)
    colors = :YlGn # [:darkred, :lightcoral, :white, :lightgreen, :forestgreen]
    return get_colorbar_normalized(vmin, vmax; colors=colors)
end

function get_colorbar_normalized(vmin::Real, vmax::Real; colors, vmid=0, rev=false)
    buckets = [vmin, vmin/2, vmid, vmax/2, vmax] # shift colormap so 0 is at center
    normed = (buckets .- vmin) / (vmax - vmin)
    return cgrad(colors, normed; rev=rev)
end

end

# ╔═╡ 5f780fbe-30bc-4ff8-bcbc-ef603a75a54c
plot_biases(results_regret, [:blob, :ellipse, :circle])

# ╔═╡ 73e83611-a98a-4bd0-a463-b0daa7529c76
plot_variances(results_regret, [:blob, :ellipse, :circle])

# ╔═╡ ec019618-4e17-4e1a-a859-c976b01a20f2
plot_mse(results_regret, [:blob, :ellipse, :circle])

# ╔═╡ 47361765-238a-414d-ab0f-ca0afb31b0cb
plot_bores(results_regret, [:blob, :ellipse, :circle])

# ╔═╡ 9f68eb56-5caf-4bf2-b015-b1d15f58109e
plot_returns(results_regret, [:blob, :ellipse, :circle])

# ╔═╡ 383c56f8-fb5d-4570-8e53-b9dae53dee33
plot_returns_var(results_regret, [:blob, :ellipse, :circle])

# ╔═╡ Cell order:
# ╟─0239ffc0-3b91-4441-86ae-cf505742d810
# ╠═e4de635f-3564-4044-a62d-d7ee9a6ba95c
# ╠═1b6df2b6-6536-4bec-aa78-39dfb90c7ad6
# ╟─a67c7a80-8db6-4d6d-b819-1ee46ebc0396
# ╠═5aa567b3-04fa-4aa5-a690-4fb93ccda580
# ╠═095f3331-54e7-472f-aea9-ab32fe45e680
# ╟─602959ed-946d-48b5-b7d8-ef56e644441c
# ╠═ec1e38f4-6e90-4041-8a96-b1d829de193c
# ╠═9085cf37-0390-482f-94b4-40e46ce3d51e
# ╠═13fc2600-3789-4fc7-a31e-357f35db4b37
# ╟─3085e225-59b9-4f0c-9586-5d82f307bf34
# ╠═b271a899-2831-456d-b5ad-5868a5890bab
# ╠═f0f3efb8-e52d-4552-9675-762c4f3fedc4
# ╠═1192839e-af4e-4ba3-9e3b-a6ebc2fedff6
# ╟─b8a2f550-43f8-4379-8cec-a0c4522c2ab9
# ╟─4e4f6167-875a-41a7-b84a-53d65f5c00a7
# ╠═620df0d7-3034-4d90-b199-84f323903f22
# ╠═b111b064-56a0-42fe-8379-155dadd80f04
# ╟─e00d2a2e-2258-4121-b96d-a80e4f30d7e1
# ╠═93545762-a44d-41a8-8e00-e1493d0748e4
# ╠═378681c1-8754-4381-bf94-fa956dc58e70
# ╟─45ad5118-df2e-44b9-9a24-20067944b06c
# ╟─ce02b5d0-4f40-4dcc-897d-8f49082725af
# ╠═80afea6c-68b1-4f4d-96c6-bae16f3bc8ac
# ╠═616e0bb0-3944-4f33-9216-afd40015382d
# ╠═7229656e-3525-4be4-a11f-dd277c3b46cd
# ╠═02c3f953-69b7-4419-96d2-18341d8f47f8
# ╠═f2be7470-82d6-4133-821a-68c887a4a2f0
# ╠═60097607-f34d-4a95-a843-efeac51297a3
# ╠═eccd50b2-4401-4cde-b2b6-021ca5263878
# ╠═0326ec74-2f6d-4d9b-93e1-fa5a0f3422a8
# ╠═29e760d4-3e7d-49ca-9698-d9c12619f720
# ╠═91c09551-3a27-4fcc-b5d7-811a652fc7db
# ╠═54cd37cf-9e53-43f5-bc3e-237cd7942c00
# ╟─e9854a55-31ad-4d7f-8331-e5575a6375c6
# ╟─ce739257-3473-426f-aef4-86a7e016a575
# ╠═44df4b43-e782-47fa-8c64-f48162aa4bd8
# ╠═cf562730-8ac6-4b45-a311-a8208c3982fb
# ╟─83a79501-d3ba-4adf-b183-4d1b28b815bd
# ╠═a3ba5f18-7d8e-4a57-85e1-56933771bd77
# ╠═9e716e52-b402-442d-9ee4-36612d2d72d0
# ╟─53d2a31e-0fd5-49a9-ae81-5a55b22e328d
# ╠═157c59f3-8034-4575-9e96-3ef3025170be
# ╠═eaebf4df-c915-4136-aa99-7101563d4aca
# ╠═74fdc433-91ff-4c08-9c46-2563be98f317
# ╠═b3296ce7-f606-4c43-b10a-3603c887be4e
# ╠═94b106f6-409d-457f-a4a7-4ad68a9fd0ef
# ╠═3bad29de-edc2-4248-b5fb-5c1aa360c457
# ╠═c1bccfa7-9f45-4ed2-8a4b-5a269329ed0a
# ╠═9b3036a7-26ac-459f-9150-3cf537258744
# ╟─c21fe459-c5c6-4324-ac17-005dddc94bb7
# ╠═e92ed5d4-955f-4138-88df-878448f1bd72
# ╠═df4bdef7-6a67-4933-b2a5-121feff6c0ee
# ╠═0a2c5707-6b04-48f9-ac2d-200392a62752
# ╠═d76fd9c4-bc40-47f9-b6c5-7334e5d2a997
# ╠═cd15ed2f-a280-4bf3-9d7e-3bb9afa4720d
# ╠═49cfde98-23b4-48d0-bed3-d739ad3b79f2
# ╟─1d2062bc-63ee-4103-b476-e1b7aab6c453
# ╠═bc1948a9-b5f2-429d-8078-c1baca0e3482
# ╠═98c42154-51c7-4274-bd46-3d53166f0575
# ╠═34bec962-d352-4bf1-a73d-993943def892
# ╠═7344ad53-3823-4a09-82b3-2600b78f66d8
# ╠═0c48c759-9203-47ee-8442-74ab7b3d5746
# ╠═1b92eaa2-63c8-4cb9-9819-05284eda89fd
# ╠═8405e36f-b6cf-481c-bb51-b67b67372c92
# ╠═41df3571-92bb-4520-b021-73dfd695d07d
# ╠═c8c7b01a-a43c-49a9-b2fd-6507b97d4130
# ╠═44dda21d-0fdd-438e-aefc-52d25397b3f2
# ╟─82bbf51b-347e-4b8d-b713-c87779eb7b3d
# ╠═97cfc22a-e2fc-4ea4-b8a5-ffbeb3f4615a
# ╠═3108dcb0-c1cb-4597-8e23-8652fede05a7
# ╠═1a5cd3f9-ad73-49e3-baf2-cf82e1a5726f
# ╠═c2110562-836f-4e4d-a999-b0fffe211da0
# ╠═8813ffa1-9a48-47c1-a096-5f88dd9b4a7c
# ╠═0ae34129-c939-4d2a-8e52-519c806ff5a5
# ╠═d3542175-9035-4dba-b2dd-257848e07510
# ╠═62e22eea-a560-4ba0-a69a-cea65301e669
# ╠═c66a8213-6d50-4bfa-b9d3-6e0e03b8b948
# ╠═e784ecf7-d5b4-4b57-93db-db019f69671a
# ╠═932679df-b6e9-416d-89bb-f9e130e6a1f7
# ╠═c95a26dd-2600-4a0e-acf2-9804cac84e8b
# ╠═fce48123-21bc-4312-a76b-86592d436dc4
# ╠═34e518bb-6c7c-4dff-ab22-7356dfd5f6fb
# ╠═61252668-bac4-4fdc-9cef-6f00fe7cb235
# ╠═f09b7aa9-d48a-465c-a246-094f168283fd
# ╠═a3c97203-6bf6-4903-ae26-eb2ed326ed08
# ╠═55ef731c-fab2-4371-9eea-f6df5cb16393
# ╠═75678d11-d516-48b5-af1b-2b31e4c5ce83
# ╟─81eb626e-3e37-47d1-b64e-8cb9e5adc9fd
# ╠═09533724-55ea-46a5-b016-3b558d452c07
# ╠═a924217b-127a-47e0-9f4e-af4e38f40ecb
# ╠═f60da90b-daf6-49ff-991a-4e7233d3b56b
# ╠═fea21277-d833-4828-b980-70f550776a01
# ╠═ec47b44d-5d53-4081-8cf7-a9cf0e923c88
# ╠═749ced45-e465-407f-b4c1-733e03802e6b
# ╠═ca17880a-ceca-42a4-82e4-6ccc16590278
# ╠═f97c751a-5ebd-4233-85a6-5a7dc59dcd09
# ╠═189c3d48-7f25-43d6-b10d-b51493b8af7e
# ╠═65d715a5-2124-4388-b901-8419a0e786a4
# ╠═33db8806-1e0c-4257-80f8-a4e4a5f2cff5
# ╠═5232cdc7-5947-4dd7-a054-955f303983ca
# ╠═14a7c259-f243-4242-8913-34deba5c0659
# ╠═eb73d59a-631d-4bca-8f94-a47f6f6137ef
# ╠═04cac42a-a4f9-4c0f-b218-bfe199c4a955
# ╠═5f780fbe-30bc-4ff8-bcbc-ef603a75a54c
# ╠═73e83611-a98a-4bd0-a463-b0daa7529c76
# ╠═ec019618-4e17-4e1a-a859-c976b01a20f2
# ╠═47361765-238a-414d-ab0f-ca0afb31b0cb
# ╠═9f68eb56-5caf-4bf2-b015-b1d15f58109e
# ╠═383c56f8-fb5d-4570-8e53-b9dae53dee33
# ╠═0cd25377-6711-47b3-9c33-35d7df0f0640
