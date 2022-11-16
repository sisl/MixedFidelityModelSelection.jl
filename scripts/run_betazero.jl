include("parallel.jl")

n = 10_000
num_cpus = 63
results_dir = abspath("MEParallel.jl/results/betazero")
config_fn = ()->sample_training_configurations(n)
results = kickoff(config_fn, num_cpus; results_dir)
