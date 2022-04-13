function goal_programming(𝐱::Dict, 𝐲::Dict, goal=[0,0]; p=2)
    return argmin(kv->norm(last(kv) - goal, p), merge(vcat, 𝐱, 𝐲))
end


function plot_pareto(results;
        minutes=true, show_origin=false, return_optimal=false, p=Inf, fn=mean)
    plot()
    mean_timings = Dict()
    std_timings = Dict()
    mean_regrets = Dict()
    std_regrets = Dict()
    std_scale = 10
    for res in results
        key, value = first(res), last(res)
        shape = key[1]
        xy = key[2][1]
        iters = key[3]

        time_divisor = minutes ? 60 : 1
        times = map(t->t.time/time_divisor, value[:timing])
        mean_time = mean(times)
        std_time = std(times)

        regrets = regret(value)
        mean_regret = fn(regrets)
        std_regret = std(regrets)

        key = "$shape, $(xy)x$(xy), $iters"
        mean_timings[key] = mean_time
        mean_regrets[key] = mean_regret
        std_timings[key] = std_time
        std_regrets[key] = std_regret

        if shape == :blob
            mark = :star6
        elseif shape == :ellipse
            mark = :diamond
        elseif shape == :circle
            mark = :circle
        end
        if xy == 10
            c = :red
        elseif xy == 30
            c = :green
        elseif xy == 50
            c = :blue
        end
        if iters == 100
            ms = 4
        elseif iters == 1000
            ms = 8
        elseif iters == 10000
            ms = 16
        end

        label = iters == 100 && shape == :circle ? "$(xy)x$(xy)" : false
        scatter!([mean_time], [mean_regret],
                 xerr=std_time/std_scale, yerr=std_regret/std_scale, msc=c,
                 c=c, ms=ms, alpha=0.3, label=label, marker=mark)
    end
    # mean_timings, mean_regrets
    pareto_optimal = goal_programming(mean_timings, mean_regrets; p=p)
    pareto_optimal_value = last(pareto_optimal)
    plot!([0, pareto_optimal_value[1]], [0, pareto_optimal_value[2]],
           c=:black, lw=1, style=:dash, label="pareto optimal")
    plot!(legend=:topright)
    title!("pareto curve (error σ/$std_scale)")
    xlabel!("mean runtime ($(minutes ? "mins." : "secs."))")
    ylabel!("mean regret")
    xlims!(0,xlims()[end])
    ylims!(0,ylims()[end])
    if show_origin
        hline!([0], c=:gray, label=false)
        vline!([0], c=:gray, label=false)
    end
    if return_optimal
        return (plot!(), pareto_optimal)
    else
        return plot!()
    end
end
