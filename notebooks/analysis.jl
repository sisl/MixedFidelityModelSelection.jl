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
results_regret = BSON.load("..\\scripts\\MEParallel.jl\\results\\results_regret.bson")[:results]

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
plot_sweep_regret(results_regret, shapekeys, regret_fn_mean, regret_title_mean)

# ╔═╡ 94b106f6-409d-457f-a4a7-4ad68a9fd0ef
begin
	α_quantile = 0.8
	regret_fn_quant = res->quantile(regret(res), α_quantile)
	regret_title_quant = "$(round(Int, 100*α_quantile))th percentile of regret"
	p_regret80 = plot_sweep_regret(results_regret, shapekeys, regret_fn_quant, 
		regret_title_quant; colorbar_regret_fn=regret_fn_90th_quant)
end

# ╔═╡ 157c59f3-8034-4575-9e96-3ef3025170be
plot_cumulative_regret(results_regret, (:blob, (10,10,1), 10); α_quantile=α_quantile)

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
p_pareto, pareto_optimal = plot_pareto(results_regret; minutes=false, return_optimal=true);

# ╔═╡ df4bdef7-6a67-4933-b2a5-121feff6c0ee
p_pareto

# ╔═╡ 0a2c5707-6b04-48f9-ac2d-200392a62752
pareto_optimal

# ╔═╡ 1d2062bc-63ee-4103-b476-e1b7aab6c453
md"""
# Solution ideas
"""

# ╔═╡ 98c42154-51c7-4274-bd46-3d53166f0575


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
# ╟─1d2062bc-63ee-4103-b476-e1b7aab6c453
# ╠═98c42154-51c7-4274-bd46-3d53166f0575
