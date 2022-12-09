include("parallel.jl")

n = 1_000
num_cpus = 63
results_dir = abspath("MEParallel.jl/results/betazero_pomcpow_obs_1k")
config_fn = ()->sample_training_configurations(n; planning_iters=1000, use_mcts=false)
results = kickoff(config_fn, num_cpus; results_dir)
