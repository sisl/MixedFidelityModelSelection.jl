# NOTE: julia1.7
using Revise
using MEParallel
using POMDPModelFidelityFramework
using Statistics
using Printf
using Plots; default(fontfamily="Palatino Roman", framestyle=:box)

if !@isdefined(results)
    @info "Loading results..."
    results = MEParallel.BSON.load(joinpath(@__DIR__, "MEParallel.jl", "results", "blob_clamped", "results_blob_clamped.bson"))[:results]
end

f(res, extraction_cost=150) = res[:r_massive]
f̂(res) = f(res) .+ last.(res[:rel_errors])

bias(res) = mean(mean(f̂(res)) .- f(res))
variance(res) = mean((f̂(res) .- mean(f̂(res))).^2)
mse(res) = mean((f(res) - f̂(res)).^2)
stddev(res) = sqrt(variance(res))
bores(res) = mean(res[:n_drills])
returns(res) = mean(res[:discounted_return])
returns_var(res) = var(res[:discounted_return])
regret_mean(res) = mean(regret(res))
regret_quantile(res) = quantile(regret(res), 0.9)
runtime(res) = mean(map(t->t.time/60, res[:timing]))

plot_biases(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, bias; title="bias", kwargs...)
plot_variances(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, variance; title="variance", kwargs...)
plot_stddev(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, stddev; title="standard deviation", kwargs...)
plot_mse(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, mse; title="mean squared error (MSE)", kwargs...)
plot_bores(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, bores; title="number of actions (drills)", kwargs...)
plot_returns(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, returns; title="mean discounted return", kwargs...)
plot_returns_var(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, returns_var; title="discounted return (variance)", kwargs...)
plot_regret(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, mean_regret; title="regret", kwargs...)
plot_regret_quantile(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, regret_quantile; title="90th percentile of regret", kwargs...)
plot_runtime(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, runtime; title="expected runtime (min.)", kwargs...)
plot_accuracy(results, shapekeys; kwargs...) = plot_fidelities(results, shapekeys, accuracy; title="accuracy", kwargs...)

function plot_fidelities(results, shapekeys, value_fn; colorbar_fn=value_fn, title="bias", titlefontsize=10, kwargs...)
    shared_cmap, svmin, svmax = get_colorbar_bounded_multi(results; value_fn=colorbar_fn)
    clims = (svmin, svmax)
    levels = 15

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
    elseif value_fn == regret_mean || value_fn == regret_quantile
        shared_cmap = cgrad(:inferno)
    elseif value_fn == runtime
        shared_cmap = cgrad(:jet)
    elseif value_fn == bores
        shared_cmap = cgrad(:jet)
    elseif value_fn == accuracy
        shared_cmap = get_colorbar_accuracy()
        clims = (0,1) # important for colormap to span [0, 1] for [red, green]
        levels = 40
    end

    is_single_belief = length(shapekeys) == 1

    pbfunc = (s; kwargs...) -> begin
        plot_fidelity(results, shapekeys[s]; reduced=true, value_func=value_fn, title=is_single_belief ? "$title" : "$(shapekeys[s])",
                  cmap=shared_cmap, clims=clims, levels=levels, titlefontsize=titlefontsize, kwargs...)
    end

    ptitle = plot(grid=false, axis=false, tick=nothing, top_margin=-5Plots.mm)

    if !is_single_belief
        title!("$title")
    end

    pplot = plot([pbfunc(s; show_ylabel=(s==1), show_xlabel=(s==2) || is_single_belief)
                 for s in 1:length(shapekeys)]...,
        layout=@layout[A B C])
    
    pcbar = contourf([0], [0], (x,y)->0,
              clims=clims, levels=levels,
              c=shared_cmap, axis=false, tick=nothing, label=false)
    # pcbar = scatter([0], [0], alpha=0,
    #     zcolor=[clims...], clims=clims, c=shared_cmap,
    #     axis=false, tick=nothing, label=false)

    layout = @layout([a{0.01h}; b c{0.1w}])

    plt_size = (735,250)
    if is_single_belief
        bottom_margin = 1Plots.mm
        left_margin = 1Plots.mm
        plt_size = (1.2*plt_size[1]÷3, 1.1*plt_size[2])
    else
        bottom_margin = 6Plots.mm
        left_margin = 3Plots.mm        
    end

    plot(ptitle, pplot, pcbar,
        formatter=format_scientific_notation,
         layout=layout,
         titlefontsize=titlefontsize,
         size=plt_size,
         bottom_margin=bottom_margin,
         left_margin=left_margin)
end


function format_scientific_notation(num)
    str = @sprintf("%.1E", num)
    split = match(r"(\d\.\d)E\+(\d+)", str)
    tens = Int(parse(Float64, split.captures[1]) * 10)
    exponent = parse(Int, split.captures[2])
    if exponent == 1
        try
            return Int(num)
        catch err
            if err isa InexactError
                return string(num)
            else
                throw(err)
            end
        end
    else
        return "$tens^{ $exponent}" # Note: en quad space " "
    end
end


function plot_fidelity(results, shapekey;
        reduced=false, cmap=nothing, clims=nothing, levels=15, title="$shapekey", titlefontsize=10,
        show_cbar=!reduced, show_xlabel=!reduced, show_ylabel=!reduced,
        value_func=res->bias(res), tight_limits=false)
    plot()
    if reduced
        title!(title)
    else
        title!("? [$title]")
    end

    if show_xlabel
        xlabel!("grid dimensions (n×n)")
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
              c=cmap, cbar=show_cbar, clims=clims, levels=levels,
              size=(500,400), yaxis=:log, linewidth=0.25, linecolor=:black)

    saved_lims = (xlims(), ylims())
    G = [(x, y) for x in X, y in Y]
    _X, _Y = first.(G), last.(G)
    scatter!(_X, _Y, ms=4, color=:white, label=false, marker=:circle)


    if tight_limits
        xlims!(saved_lims[1]...)
        ylims!(saved_lims[2]...)
    end

    return plot!(
        titlefontsize=titlefontsize,
        labelfontsize=titlefontsize-2,
        framestyle=:box,
        grid=false,
        x_foreground_color_border=:white,
        y_foreground_color_border=:white,
        x_foreground_color_axis=:white, # tick marks
        y_foreground_color_axis=:white,
    )
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
    colors = :jet # :viridis
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


function regenerate_all_heatmaps(results)
    belief_shapes = [:circle, :ellipse, :blob] # Notice low-to-high ordering
    
    plot_returns(results, belief_shapes)
    bettersavefig("heatmap_returns")
    
    plot_regret_quantile(results, belief_shapes)
    bettersavefig("heatmap_regret")
    
    plot_runtime(results, belief_shapes)
    bettersavefig("heatmap_runtime")
    
    plot_biases(results, belief_shapes)
    bettersavefig("heatmap_bias")

    plot_bores(results, belief_shapes)
    bettersavefig("heatmap_num_actions")

    plot_accuracy(results, belief_shapes)
    bettersavefig("heatmap_accuracy")

    plot_runtime(results, [:blob])
    bettersavefig("heatmap_example")
end
