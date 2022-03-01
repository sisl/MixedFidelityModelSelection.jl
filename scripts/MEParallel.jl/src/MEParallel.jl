module MEParallel

using Revise
using Distributed
using POMDPs
using POMCPOW
using Plots; default(fontfamily="Computer Modern", framestyle=:box)
using Parameters
using Statistics
using Random
using Distributions
using MineralExploration
using MixedFidelityModelSelection
using ProgressMeter
using BSON
import BSON: @save

include("interface_mineral_exploration.jl")
include("parallel_mineral_exploration.jl")

export
    # "interface_mineral_exploration.jl"
    MEJobParameters,
    MEConfiguration,
    METrial,
    MEResults,
    convert,
    job,
    mainbody_type_string,
    fileformat,
    results_filename,
    # "parallel_mineral_exploration.jl"
    configurations_fixed_bores,
    configurations_regret,
    configurations_regret20,
    makebatches,
    kickoff,
    reduce_results

end # module
