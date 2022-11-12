function regret(results; extraction_cost=150)
    massive_ore = results[:r_massive]
    r_best = max.(0, massive_ore .- extraction_cost)
    returns = results[:discounted_return]
    return r_best .- returns
end


function get_colorbar_regret(vmin, vmax)
    return cgrad(:jet)
end


function plot_sweep_regret(results, shapekeys, regret_fn, regret_title; colorbar_regret_fn=regret_fn, kwargs...)
    shared_cmap_regret, svmin_reget, svmax_regret = get_colorbar_bounded_multi(results; value_fn=colorbar_regret_fn, lower=0)

    regret_clims = (svmin_reget,svmax_regret)
    rgfunc = (s; kwargs...) -> begin
        plot_regret(results, shapekeys[s]; reduced=true, value_func=regret_fn,
            cmap=shared_cmap_regret, clims=regret_clims, kwargs...)
    end

    rgtitle = plot(title="$regret_title",
                   grid=false, axis=false, tick=nothing, top_margin=-5Plots.mm)

    rgplot = plot([rgfunc(s; show_ylabel=(s==1), show_xlabel=(s==2))
                   for s in 1:length(shapekeys)]...,
        layout=@layout[A B C])
    rgcbar = scatter([0], [0], alpha=0,
        zcolor=[regret_clims...], clims=regret_clims, c=shared_cmap_regret,
        axis=false, tick=nothing, label=false)
    rgcbar = contourf([0], [0], (x,y)->0,
        clims=regret_clims, levels=15,
        c=shared_cmap_regret, axis=false, tick=nothing, label=false)
    # rgcbar = scatter([0], [0], alpha=0,
    #     zcolor=[regret_clims...], clims=regret_clims, c=shared_cmap_regret,
    #     axis=false, tick=nothing, label=false)
    plot(rgtitle, rgplot, rgcbar,
         layout=@layout([a{0.01h}; b c{0.1w}]),
         size=(710,250), bottom_margin=6mm, left_margin=2mm)
end


function plot_regret(results, shapekey;
        reduced=false, cmap=nothing, clims=nothing,
        is_relative=false, shapekey_relative=:circle,
        show_cbar=!reduced, show_xlabel=!reduced, show_ylabel=!reduced,
        value_func=(res->mean(regret(res))), kwargs...)
    plot()
    title = "$shapekey"
    if reduced
        title!(title)
    else
        title!("expected regret [$title]")
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
        cmap = get_colorbar_regret(nothing, nothing)
    end

    if is_relative
        Z = (x,y)->abs(value_func(results[(shapekey, (x,x,1), y)]) - value_func(results[(shapekey_relative, (x,x,1), y)]))
    else
        Z = (x,y)->value_func(results[(shapekey, (x,x,1), y)])
    end

    contourf!(X, Y, Z;
              c=cmap, cbar=show_cbar, clims=clims, levels=15,
              size=(500,400), yaxis=:log, linewidth=0.25, linecolor=:black, kwargs...)

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


function plot_cumulative_regret(results, key; α_quantile=0.8)
    regrets = regret(results[key])
    shape = key[1]
    xy = key[2][1]
    iters = key[3]
    plot(x->ecdf(regrets)(x), label=false, fill=true, 0, sum(regrets), size=(300,200), fillcolor=:gray, c=:black, title="cumulative regret ($shape, $xy×$xy, $iters iters.)", titlefontsize=10)
    cr_xlims = xlims()
    cr_ylims = ylims()
    low_quant_regret = quantile(regrets, α_quantile)
    plot!([cr_xlims[1], low_quant_regret], [α_quantile, α_quantile], lw=4, c=:crimson, lab=false)
    plot!([low_quant_regret,low_quant_regret], [0, α_quantile], lw=4, c=:crimson, lab=false)
    scatter!([low_quant_regret], [α_quantile], ms=4, c=:red, lab=false)
    xlims!(cr_xlims...)
    ylims!(0, cr_ylims[2])
end
