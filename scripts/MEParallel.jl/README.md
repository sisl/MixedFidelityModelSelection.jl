# MEParallel.jl

Launch julia with `julia -p auto`.

```julia
include("parallel.jl")

kickoff(()->configurations_standard("inject_perturbed"), 63; results_dir=abspath("MEParallel.jl/results/inject_perturbed"))
```