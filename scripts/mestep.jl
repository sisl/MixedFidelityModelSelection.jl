using Revise
using Infiltrator
using Plots; default(fontfamily="Computer Modern", frame=:box)
using POMDPSimulators

using MineralExploration
include("MEParallel.jl")
using .MEParallel

config = MEConfiguration(1, (30,30,1), 10, CircleNode, false, MEJobParameters(name="step", min_bores=5, max_bores=5))
trial = initialize(config)

display(plot(trial.b0))
for (sp,a,o,bp,t) in stepthrough(trial.problem, trial.planner, trial.updater, trial.b0, trial.s0, "sp,a,o,bp,t", max_steps=5)
    @info "Time $t"
    display(plot(bp))
    @infiltrate
end
