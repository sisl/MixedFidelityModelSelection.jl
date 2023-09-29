# POMDPModelFidelityFramework.jl

Model-fidelity analysis for POMDPs, with a case study for critical mineral exploration ([MineralExploration.jl](https://github.com/sisl/MineralExploration)).

<p align="center">
    <img src="./media/PMFF.jpg">
</p>

## Interface

```julia
# (::Type{<:Configuration})::Vector{<:Configuration}
function configurations end

# (config::Configuration)::Trial
function initialize end

# (trial::Trial)::Results
function evaluate end

# (results::Results)
function save end
```

To see an example of an implemented interface, see: [`interface_mineral_exploration.jl`](https://github.com/sisl/POMDPModelFidelityFramework.jl/blob/main/scripts/MEParallel.jl/src/interface_mineral_exploration.jl)

---

Contact: [Robert Moss](https://github.com/mossr)

