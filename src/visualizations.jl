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
