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