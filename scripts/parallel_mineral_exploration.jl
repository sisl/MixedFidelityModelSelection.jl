using Distributed

using Revise
using POMDPs
using POMCPOW
using Plots; default(fontfamily="Computer Modern", framestyle=:box)
using Parameters
using Statistics
using Random
using Distributions
using MineralExploration
using MixedFidelityModelSelection

@with_kw struct MEJobParameters <: JobParameters
    high_fidelity_dim::Tuple{Real,Real,Real} = (50,50,1)
    min_bores::Int = 5
    max_bores::Int = 25
    grid_spacing::Int = 0
    max_movement::Int = 20
    n_initial::Int = 0
end

struct MEConfiguration <: Configuration
    seed::Int
    grid_dims::Tuple{Real,Real,Real}
    pomcpow_iters::Int
    mainbody_type::Type{<:MainbodyGen}
    params::MEJobParameters
end

struct METrial <: Trial
    config::MEConfiguration
    problem::Union{MDP, POMDP}
    planner::Policy
    s0 # initial state
    b0 # initial belief
    updater::Updater
end

@with_kw struct MEResults <: Results 
    config::MEConfiguration
    timing::NamedTuple
    ore_map::Array
    massive_threshold::Real

    # from `run_trial` in utils.jl
    discounted_return
    distances
    abs_errors
    rel_errors
    vol_stds
    n_drills
    r_massive
    last_action
end

function configurations(::Type{MEConfiguration})
    params = MEJobParameters(min_bores=20, max_bores=20)
    configs = MEConfiguration[]
    for seed in 1:10
        for pomcpow_iters in [100]
            for grid_dims_xy in [30, 50]
                for mainbody_type in [BlobNode, EllipseNode, CircleNode]
                    grid_dims = (grid_dims_xy, grid_dims_xy, 1)
                    config = MEConfiguration(seed, grid_dims, pomcpow_iters, mainbody_type, params)
                    push!(configs, config)
                end
            end
        end
    end
    return configs::Vector{<:Configuration}
end

function initialize(config::Configuration)
    # Configuration specific option
    seed = config.seed
    grid_dims = config.grid_dims
    mainbody = config.mainbody_type(grid_dims=grid_dims)
    pomcpow_iters = config.pomcpow_iters

    # Static job parameters
    high_fidelity_dim = config.params.high_fidelity_dim
    min_bores = config.params.min_bores
    max_bores = config.params.max_bores
    grid_spacing = config.params.grid_spacing
    max_movement = config.params.max_movement
    n_initial = config.params.n_initial

    is10x10 = grid_dims == (10,10,1)
    original_max_movement = is10x10 ? 0 : max_movement
    grid_spacing_delta = is10x10 ? 0 : grid_spacing+1

    Random.seed!(seed)
    m = MineralExplorationPOMDP(grid_dim=grid_dims,
                                high_fidelity_dim=high_fidelity_dim,
                                delta=grid_spacing_delta,
                                grid_spacing=grid_spacing,
                                mainbody_gen=mainbody,
                                original_max_movement=original_max_movement,
                                min_bores=min_bores,
                                max_bores=max_bores)
    initialize_data!(m, n_initial)

    Random.seed!(seed)
    ds0 = POMDPs.initialstate_distribution(m)
    s0 = rand(ds0; truth=true)

    Random.seed!(seed)
    up = MEBeliefUpdater(m, 1000, 2.0)
    b0 = POMDPs.initialize_belief(up, ds0)

    solver = POMCPOWSolver(tree_queries=pomcpow_iters,
                           check_repeat_obs=true,
                           check_repeat_act=true,
                           next_action=NextActionSampler(),
                           k_action=2.0,
                           alpha_action=0.25,
                           k_observation=2.0,
                           alpha_observation=0.1,
                           criterion=POMCPOW.MaxUCB(100.0),
                           final_criterion=POMCPOW.MaxQ(),
                           # final_criterion=POMCPOW.MaxTries(),
                           estimate_value=0.0
                           # estimate_value=leaf_estimation
                           )
    planner = POMDPs.solve(solver, m)
    trial = METrial(config, m, planner, s0, b0, up)
    return trial
end

function evaluate(trial::Trial; save_dir=nothing)
    grid_dims = trial.config.grid_dims
    pomcpow_iters = trial.config.pomcpow_iters
    seed = trial.config.seed

    m = trial.problem
    up = trial.updater
    planner = trial.planner
    s0 = trial.s0
    b0 = trial.b0

    !isnothing(save_dir) && !isdir(save_dir) && mkdir(save_dir)

    Random.seed!(trial.config.seed)
    timing = @timed trial_results = run_trial(m, up, planner, s0, b0, save_dir=save_dir, display_figs=false)
    return trial_results, m, s0.ore_map, grid_dims, pomcpow_iters, seed, timing
end






# @everywhere begin
#     N_INITIAL = 0
#     MAX_BORES = 25
#     MIN_BORES = 5
#     # MAX_BORES = 20
#     # MIN_BORES = 20
#     GRID_SPACING = 0
#     MAX_MOVEMENT = 20
#     SAVE_DIR = "../results/data/sensitivity_analysis_sweep/"
#     !isdir(SAVE_DIR) && mkdir(SAVE_DIR)
# end

# @everywhere begin
#     function run_single_trial(seed, grid_dim_xy, pomcpow_iter; TrueMainbodyType=BlobNode, MainbodyType=BlobNode, save_dir=nothing)
#         grid_dims = (grid_dim_xy, grid_dim_xy, 1)
#         @show seed, grid_dims, pomcpow_iter, TrueMainbodyType, MainbodyType
#         high_fidelity_dim = (50,50,1) # TODO.
#         true_mainbody = TrueMainbodyType(grid_dims=high_fidelity_dim)
#         mainbody = MainbodyType(grid_dims=grid_dims)

#         is10x10 = grid_dims == (10,10,1)
#         original_max_movement = is10x10 ? 0 : MAX_MOVEMENT
#         grid_spacing_delta = is10x10 ? 0 : GRID_SPACING+1

#         Random.seed!(seed)
#         m = MineralExplorationPOMDP(max_bores=MAX_BORES, delta=grid_spacing_delta,
#                                     grid_spacing=GRID_SPACING,
#                                     true_mainbody_gen=true_mainbody,
#                                     mainbody_gen=mainbody,
#                                     original_max_movement=MAX_MOVEMENT,
#                                     min_bores=MIN_BORES, grid_dim=grid_dims,
#                                     high_fidelity_dim=high_fidelity_dim)
#         initialize_data!(m, N_INITIAL)

#         Random.seed!(seed)
#         ds0 = POMDPs.initialstate_distribution(m)
#         s0 = rand(ds0; truth=true)

#         # Run experiment
#         Random.seed!(seed)
#         up = MEBeliefUpdater(m, 1000, 2.0)
#         b0 = POMDPs.initialize_belief(up, ds0)

#         next_action = NextActionSampler()
#         solver = POMCPOWSolver(tree_queries=pomcpow_iter,
#                                check_repeat_obs=true,
#                                check_repeat_act=true,
#                                next_action=next_action,
#                                k_action=2.0,
#                                alpha_action=0.25,
#                                k_observation=2.0,
#                                alpha_observation=0.1,
#                                criterion=POMCPOW.MaxUCB(100.0),
#                                final_criterion=POMCPOW.MaxQ(),
#                                # final_criterion=POMCPOW.MaxTries(),
#                                estimate_value=0.0
#                                # estimate_value=leaf_estimation
#                                )
#         planner = POMDPs.solve(solver, m)
#         timing = @timed results = run_trial(m, up, planner, s0, b0, save_dir=save_dir, display_figs=false)
#         return results, m, s0.ore_map, grid_dim_xy, pomcpow_iter, seed, timing
#     end
# end

# target_procs = 9 # TODO: better parallelization over configs and seeds (9*5)
# if nprocs() < target_procs
#     addprocs(target_procs - nprocs())
# end
# @show nprocs()

# include("run_sensitivity_analysis_sweep_parallel.jl")