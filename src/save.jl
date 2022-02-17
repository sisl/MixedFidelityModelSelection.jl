import BSON: @save

m = poutput[1][2]
ore_maps = [p[3].ore_map for p in poutput]

analysis = Dict(
    "results" => first.(poutput),
    "timing" => last.(poutput),
    "grid_dims" => [p[4] for p in poutput],
    "pomcpow_iters" => [p[5] for p in poutput],
    "seeds" => [p[6] for p in poutput],
    "ore_maps" => ore_maps,
    "mass_maps" => [ore_map .>= m.massive_threshold for ore_map in ore_maps],
    "r_massives" => [1*sum(ore_map .>= m.massive_threshold) for ore_map in ore_maps], # NOTE: dim_scale == 1
)

@save filename analysis
