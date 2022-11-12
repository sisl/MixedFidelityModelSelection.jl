include("parallel.jl")

p = SimDecParameters()
n = 500*9*3
num_cpus = 39
results_dir = abspath("MEParallel.jl/results/simdec")
results = kickoff(()->sample_simdec_configurations(p, n), num_cpus; results_dir)
save_simdec_csv(results, results_dir)
