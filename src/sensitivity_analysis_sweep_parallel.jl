using Distributed

# @everywhere using Revise
@everywhere using POMDPs
@everywhere using POMCPOW
@everywhere using Plots; default(fontfamily="Computer Modern", framestyle=:box)
@everywhere using Statistics
@everywhere using Random
@everywhere using Distributions
@everywhere using MineralExploration

@everywhere plot_ore_map = MineralExploration.plot_ore_map
@everywhere plot_mass_map = MineralExploration.plot_mass_map

@everywhere begin
    N_INITIAL = 0
    MAX_BORES = 25
    MIN_BORES = 5
    # MAX_BORES = 20
    # MIN_BORES = 20
    GRID_SPACING = 0
    MAX_MOVEMENT = 20
    SAVE_DIR = "../results/data/sensitivity_analysis_sweep/"
    !isdir(SAVE_DIR) && mkdir(SAVE_DIR)
end

@everywhere begin
    function sample_massive(m,
            ore_bounds=[m.extraction_cost-10, m.extraction_cost+10];
            N=1, max_iters=1000)
        ds0 = POMDPs.initialstate_distribution(m)
        samples = []
        seeds = []
        for _ in 1:max_iters
            seed = rand(0:typemax(Int32))
            Random.seed!(seed)
            s0 = rand(ds0; truth=true)
            s_massive = s0.ore_map .>= m.massive_threshold
            r_massive = m.dim_scale*sum(s_massive)
            if ore_bounds[1] ≤ round(r_massive, digits=2) ≤ ore_bounds[2]
                push!(samples, s0.ore_map)
                push!(seeds, seed)
            end
            length(samples) ≥ N && break
        end
        return samples, seeds
    end
end

@everywhere begin
    function run_single_trial(seed, grid_dim_xy, pomcpow_iter; TrueMainbodyType=BlobNode, MainbodyType=BlobNode, save_dir=nothing)
        grid_dims = (grid_dim_xy, grid_dim_xy, 1)
        @show seed, grid_dims, pomcpow_iter, TrueMainbodyType, MainbodyType
        high_fidelity_dim = (50,50,1) # TODO.
        true_mainbody = TrueMainbodyType(grid_dims=high_fidelity_dim)
        mainbody = MainbodyType(grid_dims=grid_dims)

        is10x10 = grid_dims == (10,10,1)
        original_max_movement = is10x10 ? 0 : MAX_MOVEMENT
        grid_spacing_delta = is10x10 ? 0 : GRID_SPACING+1

        Random.seed!(seed)
        m = MineralExplorationPOMDP(max_bores=MAX_BORES, delta=grid_spacing_delta,
                                    grid_spacing=GRID_SPACING,
                                    true_mainbody_gen=true_mainbody,
                                    mainbody_gen=mainbody,
                                    original_max_movement=MAX_MOVEMENT,
                                    min_bores=MIN_BORES, grid_dim=grid_dims,
                                    high_fidelity_dim=high_fidelity_dim)
        initialize_data!(m, N_INITIAL)

        Random.seed!(seed)
        ds0 = POMDPs.initialstate_distribution(m)
        s0 = rand(ds0; truth=true)

        # Run experiment
        Random.seed!(seed)
        up = MEBeliefUpdater(m, 1000, 2.0)
        b0 = POMDPs.initialize_belief(up, ds0)

        next_action = NextActionSampler()
        solver = POMCPOWSolver(tree_queries=pomcpow_iter,
                               check_repeat_obs=true,
                               check_repeat_act=true,
                               next_action=next_action,
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
        timing = @timed results = run_trial(m, up, planner, s0, b0, save_dir=save_dir, display_figs=false)
        return results, m, s0.ore_map, grid_dim_xy, pomcpow_iter, seed, timing
    end
end

target_procs = 9 # TODO: better parallelization over configs and seeds (9*5)
if nprocs() < target_procs
    addprocs(target_procs - nprocs())
end
@show nprocs()

# include("run_sensitivity_analysis_sweep_parallel.jl")