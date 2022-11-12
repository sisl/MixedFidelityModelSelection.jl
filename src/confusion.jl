count_true_positives(decisions, true_decisions) = sum(decisions .== :mine .&& true_decisions .== :mine)
count_true_negatives(decisions, true_decisions) = sum(decisions .== :abandon .&& true_decisions .== :abandon)
count_false_positives(decisions, true_decisions) = sum(decisions .== :mine .&& true_decisions .== :abandon)
count_false_negatives(decisions, true_decisions) = sum(decisions .== :abandon .&& true_decisions .== :mine)


function Base.precision(decisions::Vector, true_decisions::Vector)
    tp = count_true_positives(decisions, true_decisions)
    fp = count_false_positives(decisions, true_decisions)
    return tp / (tp+fp)
end


function recall(decisions, true_decisions)
    tp = count_true_positives(decisions, true_decisions)
    fn = count_false_negatives(decisions, true_decisions)
    return tp / (tp+fn)
end


function get_true_decisions(results::Dict; extraction_cost=150)
    return map(r_massive->r_massive > extraction_cost ? :mine : :abandon, results[:r_massive])
end


function accuracy(results::Dict; extraction_cost=150)
    decisions = results[:last_action]
    true_decisions = get_true_decisions(results; extraction_cost=extraction_cost)
    return accuracy(decisions, true_decisions)
end


function accuracy(decisions, true_decisions)
    tp = count_true_positives(decisions, true_decisions)
    tn = count_true_negatives(decisions, true_decisions)
    fp = count_false_positives(decisions, true_decisions)
    fn = count_false_negatives(decisions, true_decisions)
    return (tp+tn) / (tp+tn+fp+fn)
end


function confusion_matrix(decisions, true_decisions)
    N = length(decisions)
    tpr = count_true_positives(decisions, true_decisions) / N
    tnr = count_true_negatives(decisions, true_decisions) / N
    fpr = count_false_positives(decisions, true_decisions) / N
    fnr = count_false_negatives(decisions, true_decisions) / N
    return [tpr fpr;
            fnr tnr]
end


function confusion_table(cm, shape, grid_dims)
    title = "\$\\text{$shape: }$grid_dims\\times$grid_dims\$"
    label = " \${}^\\text{true}/_\\text{false}\$" # Note leading space (Pluto needs this).
    return Markdown.parse("""
$title
$label   | positive   | negative
:------- | :--------- | :-------
positive | $(cm[1,1]) | $(cm[1,2])
negative | $(cm[2,1]) | $(cm[2,2])
""")
end


function plot_confusion(cm)
    heatmap(rotr90(cm)', c=:viridis, ratio=1, xlims=(0.5, 2.5), ylims=(0.5, 2.5), clims=(0,1))
end


function get_colorbar_accuracy()
    return get_colorbar(0, 1; vmid=0.5)
end


function plot_accuracies(results, shapekeys)
    mbfunc = (s; kwargs...) -> begin
        plot_accuracy(results, shapekeys[s]; reduced=true, kwargs...)
    end

    mbtitle = plot(title="accuracy",
                   grid=false, axis=false, tick=nothing, top_margin=-5Plots.mm)

    mbplot = plot([mbfunc(s; show_ylabel=(s==1), show_xlabel=(s==2))
                   for s in 1:length(shapekeys)]...,
        layout=@layout[A B C], size=(700,210)
    )
    pcbar = contourf([0], [0], (x,y)->0,
              clims=(0,1), levels=15,
              c=get_colorbar(0, 1; vmid=0.5), axis=false, tick=nothing, label=false)
    # pcbar = scatter([0], [0], alpha=0,
    #     zcolor=[0,1], clims=(0,1), c=get_colorbar(0, 1; vmid=0.5),
    #     axis=false, tick=nothing, label=false)
    plot(mbtitle, mbplot, pcbar,
         layout=@layout([a{0.01h}; b c{0.1w}]),
         size=(710,250), bottom_margin=6mm, left_margin=2mm)
end


function plot_accuracy(results, shapekey;
        reduced=false, cmap=nothing,
        show_cbar=!reduced, show_xlabel=!reduced, show_ylabel=!reduced,
        value_fn=accuracy, # betavalue,
        title="accuracy",
    )
    plot()
    subtitle = "$shapekey"
    if reduced
        title!(subtitle)
    else
        title!("$title [$subtitle]")
    end

    if show_xlabel
        xlabel!("grid dimension (n×n)")
    end
    if show_ylabel
        ylabel!("planning iterations")
    end

    X, Y = get_data_xy(results)

    # Color scheme (centered at 0.5)
    if isnothing(cmap)
        cmap = get_colorbar_accuracy()
    end

    contourf!(X, Y, (x,y)->value_fn(results[(shapekey, (x,x,1), y)]),
              c=cmap, cbar=show_cbar, clims=(0,1), levels=16,
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

