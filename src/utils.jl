function sample_massive(m,
        ore_bounds=[m.extraction_cost-10, m.extraction_cost+10];
        N=1, max_iters=1000)
    ds0 = POMDPs.initialstate_distribution(m)
    samples = []
    seeds = []
    for _ in 1:max_iters
        seed = rand(0:typemax(Int32))
        Random.seed!(seed)
        s0 = rand(ds0; truth=true)
        s_massive = s0.ore_map .>= m.massive_threshold
        r_massive = m.dim_scale*sum(s_massive)
        if ore_bounds[1] ≤ r_massive ≤ ore_bounds[2]
            push!(samples, s0.ore_map)
            push!(seeds, seed)
        end
        length(samples) ≥ N && break
    end
    return samples, seeds
end


function filltomatch!(A::Vector{Vector{T}}) where T <: Any
    max_length = maximum(map(length, A))    
    for i in 1:length(A)
        len = length(A[i])
        if len < max_length
            # fill last element to `max_length`
            A[i] = vcat(A[i], fill(A[i][end], max_length-len))
        end
    end
    return A
end


function mismatch_mean(A)
    max_length = maximum(map(length, A))
    Z = [map(a->i <= length(a) ? a[i] : nothing, A) for i in 1:max_length]
    return map(mean, map(z->filter(!isnothing, z), Z))
end


function mismatch_std(A)
    max_length = maximum(map(length, A))
    Z = [map(a->i <= length(a) ? a[i] : nothing, A) for i in 1:max_length]
    stds = map(std, map(z->filter(!isnothing, z), Z))
    return map(σ->isnan(σ) ? 0 : σ, stds)
end

bettersavefig(filename; kwargs...) = bettersavefig(plot!(), filename; kwargs...)

function bettersavefig(fig, filename; dpi=300)
    filename_png = "$filename.png"
    filename_svg = "$filename.svg"
    savefig(fig, filename_svg)
    run(`inkscape -f $filename_svg -e $filename_png -d $dpi`)
    rm(filename_svg)
end
