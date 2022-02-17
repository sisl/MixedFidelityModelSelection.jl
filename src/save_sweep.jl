m = poutput[1][2][1]
ore_maps = [p[3] for p in poutput]

analysis = Dict(
    "results" => first.(poutput),
    "timing" => last.(poutput),
    "grid_dim_xys" => [p[4] for p in poutput],
    "pomcpow_iters" => [p[5] for p in poutput],
    "seeds" => [p[6] for p in poutput],
    "ore_maps" => ore_maps,
    "mass_maps" => [[ore_map[i] .>= m.massive_threshold for i in 1:length(ore_map)] for ore_map in ore_maps],
    "r_massives" => [[1*sum(ore_map[i] .>= m.massive_threshold) for i in 1:length(ore_map)] for ore_map in ore_maps], # NOTE: dim_scale == 1
)

@save filename analysis
