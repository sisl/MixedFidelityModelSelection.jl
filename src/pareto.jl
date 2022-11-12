normalize01(x, min, max) = (x - min) / (max - min)


function goal_programming(𝐱::Dict, 𝐲::Dict, goal=[0,0]; p=2)
    return argmin(kv->norm(last(kv) - goal, p), merge(vcat, 𝐱, 𝐲))
end


function goal_programming(𝐱::Dict, 𝐲::Dict, 𝐳::Dict, goal=[0,0,0]; p=2)
    return argmin(kv->norm(last(kv) - goal, p), merge(vcat, 𝐱, 𝐲, 𝐳))
end


function plot_pareto(results;
        minutes=true, show_origin=false, return_optimal=false, p=Inf, fn=mean, yfunc=regret, ylabel="mean regret", normalized=true)
    plot()
    mean_timings = Dict()
    std_timings = Dict()
    mean_yvalues = Dict()
    std_yvalues = Dict()
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

        yvalues = yfunc(value)
        mean_yvalue = fn(yvalues)
        std_yvalue = std(yvalues)

        key = "$shape, $(xy)x$(xy), $iters"
        mean_timings[key] = mean_time
        mean_yvalues[key] = mean_yvalue
        std_timings[key] = std_time
        std_yvalues[key] = std_yvalue
    end

    if normalized
        min_timings_μ, max_timings_μ = minimum(values(mean_timings)), maximum(values(mean_timings))
        min_timings_σ, max_timings_σ = minimum(values(std_timings)), maximum(values(std_timings))

        min_yvalues_μ, max_yvalues_μ = minimum(values(mean_yvalues)), maximum(values(mean_yvalues))
        min_yvalues_σ, max_yvalues_σ = minimum(values(std_yvalues)), maximum(values(std_yvalues))

        for key in keys(mean_timings)
            mean_timings[key] = normalize01(mean_timings[key], min_timings_μ, max_timings_μ)
            mean_yvalues[key] = normalize01(mean_yvalues[key], min_yvalues_μ, max_yvalues_μ)

            std_timings[key] = normalize01(std_timings[key], min_timings_σ, max_timings_σ)
            std_yvalues[key] = normalize01(std_yvalues[key], min_yvalues_σ, max_yvalues_σ)
        end 
    end 

    for res in results
        key, value = first(res), last(res)
        shape = key[1]
        xy = key[2][1]
        iters = key[3]

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

        key = "$shape, $(xy)x$(xy), $iters"
        mean_time = mean_timings[key]
        mean_yvalue = mean_yvalues[key]
        std_time = std_timings[key]
        std_yvalue = std_yvalues[key]          

        label = iters == 1000 && shape == :circle ? "$(xy)x$(xy)" : false
        scatter!([mean_time], [mean_yvalue],
                 xerr=std_time/std_scale, yerr=std_yvalue/std_scale, msc=c,
                 c=c, ms=ms, alpha=0.3, label=label, marker=mark)
    end
    # mean_timings, mean_yvalues
    pareto_optimal = goal_programming(mean_timings, mean_yvalues; p=p)
    pareto_optimal_value = last(pareto_optimal)
    plot!([0, pareto_optimal_value[1]], [0, pareto_optimal_value[2]],
           c=:black, lw=1, style=:dash, label="pareto optimal")
    scatter!([0], [0], marker=:square, label=false, color=:white)
    plot!(legend=:topright)
    title!("pareto curve (error σ/$std_scale)")
    if normalized
        xlabel!("mean runtime (normalized)")
        ylabel!("$ylabel (normalized)")
    else
        xlabel!("mean runtime ($(minutes ? "mins." : "secs."))")
        ylabel!(ylabel)
    end
    if !normalized
        xlims!(0,xlims()[end])
        ylims!(0,ylims()[end])
    end
    if show_origin
        hline!([0], c=:gray, label=false)
        vline!([0], c=:gray, label=false)
    end

    plt = normalized ? plot!(ratio=1, xlims=ylims()) : plot!()

    if return_optimal
        return (plt, pareto_optimal)
    else
        return plt
    end
end



function plot_pareto_3d(results;
        minutes=true, show_origin=false, return_optimal=false, p=Inf,
        yfn=mean, yfunc=regret, ylabel="mean regret",
        zfn=x->x, zfunc=x->1-accuracy(x), zlabel="1 - accuracy",
        normalized=true)
    plot()
    mean_timings = Dict()
    std_timings = Dict()
    mean_yvalues = Dict()
    std_yvalues = Dict()
    mean_zvalues = Dict()
    std_zvalues = Dict()
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

        yvalues = yfunc(value)
        mean_yvalue = yfn(yvalues)
        std_yvalue = std(yvalues)

        zvalues = zfunc(value)
        mean_zvalue = zfn(zvalues)
        std_zvalue = std(zvalues)

        key = "$shape, $(xy)x$(xy), $iters"
        mean_timings[key] = mean_time
        mean_yvalues[key] = mean_yvalue
        mean_zvalues[key] = mean_zvalue
        std_timings[key] = std_time
        std_yvalues[key] = std_yvalue
        std_zvalues[key] = std_zvalue
    end

    if normalized
        min_timings_μ, max_timings_μ = minimum(values(mean_timings)), maximum(values(mean_timings))
        min_timings_σ, max_timings_σ = minimum(values(std_timings)), maximum(values(std_timings))

        min_yvalues_μ, max_yvalues_μ = minimum(values(mean_yvalues)), maximum(values(mean_yvalues))
        min_yvalues_σ, max_yvalues_σ = minimum(values(std_yvalues)), maximum(values(std_yvalues))

        min_zvalues_μ, max_zvalues_μ = minimum(values(mean_zvalues)), maximum(values(mean_zvalues))
        min_zvalues_σ, max_zvalues_σ = minimum(values(std_zvalues)), maximum(values(std_zvalues))

        for key in keys(mean_timings)
            mean_timings[key] = normalize01(mean_timings[key], min_timings_μ, max_timings_μ)
            mean_yvalues[key] = normalize01(mean_yvalues[key], min_yvalues_μ, max_yvalues_μ)
            mean_zvalues[key] = normalize01(mean_zvalues[key], min_zvalues_μ, max_zvalues_μ)

            std_timings[key] = normalize01(std_timings[key], min_timings_σ, max_timings_σ)
            std_yvalues[key] = normalize01(std_yvalues[key], min_yvalues_σ, max_yvalues_σ)
            std_zvalues[key] = normalize01(std_zvalues[key], min_zvalues_σ, max_zvalues_σ)
        end 
    end 

    for res in results
        key, value = first(res), last(res)
        shape = key[1]
        xy = key[2][1]
        iters = key[3]

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

        key = "$shape, $(xy)x$(xy), $iters"
        mean_time = mean_timings[key]
        mean_yvalue = mean_yvalues[key]
        mean_zvalue = mean_zvalues[key]
        std_time = std_timings[key]
        std_yvalue = std_yvalues[key]          
        std_zvalue = std_zvalues[key]          

        label = iters == 1000 && shape == :circle ? "$(xy)x$(xy)" : false
        scatter!([mean_time], [mean_yvalue], [mean_zvalue], msc=c,
                 c=c, ms=ms, alpha=0.3, label=label, marker=mark, camera=(45,45))
    end
    # mean_timings, mean_yvalues
    pareto_optimal = goal_programming(mean_timings, mean_yvalues, mean_zvalues; p=p)
    pareto_optimal_value = last(pareto_optimal)
    plot!([0, pareto_optimal_value[1]], [0, pareto_optimal_value[2]], [0, pareto_optimal_value[3]],
           c=:black, lw=1, style=:dash, label="pareto optimal")
    plot!(legend=:outertopright)
    title!("pareto curve (error σ/$std_scale)")
    if normalized
        xlabel!("||mean runtime||")
        ylabel!("||$ylabel||")
        zlabel = "||$zlabel||"
    else
        xlabel!("mean runtime ($(minutes ? "mins." : "secs."))")
        ylabel!(ylabel)
    end
    if !normalized
        xlims!(0,xlims()[end])
        ylims!(0,ylims()[end])
    end

    plt = normalized ? plot!(ratio=1, xlims=ylims(), zlabel=zlabel) : plot!()

    if return_optimal
        return (plt, pareto_optimal)
    else
        return plt
    end
end
