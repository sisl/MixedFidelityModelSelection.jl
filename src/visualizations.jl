tufte_colors = Dict(
    "pastel_magenta" => colorant"#FF48CF",
    "pastel_purple" => colorant"#8770FE",
    "pastel_blue" => colorant"#1BA1EA",
    "pastel_seagreen" => colorant"#14B57F",
    "pastel_green" => colorant"#3EAA0D",
    "pastel_orange" => colorant"#C38D09",
    "pastel_red" => colorant"#F5615C",  
)

function plot_ore_mass_distributions(ore_masses, labels, grid_dims)
    latexize = str -> latexstring("\\mathrm{$(replace(str, " "=>"\\; "))}")

    shapes = latexize.(labels)
    
    fig = Figure(font="Computer Modern")
    offset_scale = 32
    dim_title = "$(grid_dims[1])×$(grid_dims[2])"

    # Note, reversed (highest fidelity on top)
    axes = Axis(fig[1,1], title=L"ore mass distributions $%$dim_title$",
         yticks=((1:length(shapes)) ./ offset_scale, reverse(shapes)),
         xlabel=latexize("ore mass"),
         xtickformat = xs -> [latexize(string(Int(x))) for x in xs])

    for i in 1:length(shapes)
        # Note, reversed (highest fidelity on top)
        data = convert(Vector{Real}, reverse(ore_masses)[i])
        d = CairoMakie.density!(data, offset=i/offset_scale,
                                color=:x, colormap=:solar,
                                strokewidth=1, strokecolor=:black)
        μ = mean(data)
        σ = std(data)
        Xs = [μ, μ]
        Ys = [i/offset_scale, (i+0.5)/offset_scale]
        # Ys[2] = i == 3 ? Ys[2]*1.05 : Ys[2]
        CairoMakie.lines!(Xs, Ys, color=:orange, linewidth=3)
        CairoMakie.text!(latexstring("$(round(μ, digits=2)) ± $(round(σ, digits=2))"),
                         position=(Xs[2], Ys[2]+(0.05/offset_scale)),
                         align=(:center, :baseline),
                         textsize=14)
    end
    return fig
end


draw_zero() = plot!([xlims()...], [0, 0], lw=1, c=:gray, xlims=xlims(), label=false)


function plot_relative_errors(rel_errs_vector; c=:crimson, label=false, hold=false, marker=:circle, grid_dim=NaN)
    gds = string(grid_dim, "x", grid_dim)
    μ_rel_errs = mismatch_mean(rel_errs_vector)
    σ_rel_errs = mismatch_std(rel_errs_vector)
    ts = [1:length(μ_rel_errs);] .- 1
    plot_fn = hold ? plot! : plot
    plot_fn(ts, μ_rel_errs,
            title="relative volume error (aggregate: $gds, 10 seeds)",
            ribbon=σ_rel_errs, fillalpha=0.1,
            marker=(marker, 3, 1.0),
            bg_inside=:white,
            xlabel="time step", ylabel="relative error", label=label, lw=2, c=c)
    plot!(ts, μ_rel_errs .+ σ_rel_errs, c=c, label=false, alpha=0.33)
    plot!(ts, μ_rel_errs .- σ_rel_errs, c=c, label=false, alpha=0.33)
    annotate!(0.2, 10, text("overestimated", :gray, :left, 10, "Computer Modern"))
    annotate!(0.2, -10, text("underestimated", :gray, :left, 10, "Computer Modern"))
    xlims!(0, Inf)
    # ylims!(minimum(μ_rel_errs), Inf)
    plot!(legend=:bottomright)
end


function plot_rel_error_aggregate(resultsblob, resultsellipse, resultscircle;
        mini=false)
    rel_errs_blob = resultsblob[:rel_errors]
    rel_errs_ellipse = resultsellipse[:rel_errors]
    rel_errs_circle = resultscircle[:rel_errors]
    grid_dim = resultsblob[:config][1][:grid_dims][1]
    colors = cgrad(:lighttest, 3, categorical=true)

    plot_relative_errors(rel_errs_blob;
        label=string("blob: ", accuracy(resultsblob)),
        hold=false, c=colors[1], marker=:star6, grid_dim=grid_dim)
    plot_relative_errors(rel_errs_ellipse;
        label=string("ellipse: ", accuracy(resultsellipse)),
        hold=true, c=colors[2], marker=:diamond, grid_dim=grid_dim)
    plot_relative_errors(rel_errs_circle;
        label=string("circle: ", accuracy(resultscircle)),
        hold=true, c=colors[3], marker=:circle, grid_dim=grid_dim)
    draw_zero()

    if mini
        xlims!(15, 20)
        ylims!(-35, 5)
        title!("")
        # xlabel!("")
        # ylabel!("")
        plot!(legend=false, size=(300,250))
    else
        plot!()
    end
    ylims!(-140, 140)
end


