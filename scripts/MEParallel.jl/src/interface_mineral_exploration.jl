@with_kw struct MEJobParameters <: JobParameters
    high_fidelity_dim::Tuple{Real,Real,Real} = (50,50,1)
    min_bores::Int = 5
    max_bores::Int = 25
    grid_spacing::Int = 0
    max_movement::Int = 20
    n_initial::Int = 0
    name::String = "job"
end


struct MEConfiguration <: Configuration
    seed::Int
    grid_dims::Tuple{Real,Real,Real}
    pomcpow_iters::Int
    mainbody_type::Type{<:MainbodyGen}
    use_mcts::Bool
    params::MEJobParameters
end


struct METrial <: Trial
    config::MEConfiguration
    problem::Union{MDP, POMDP}
    planner::Policy
    s0::Any # initial state
    b0::Any # initial belief
    updater::Updater
end


@with_kw struct MEResults <: Results
    config::MEConfiguration
    timing::NamedTuple
    ore_map::Array{Real,3}
    mass_map::Array{Real,3}
    massive_threshold::Real

    # from `run_trial` in utils.jl
    discounted_return::Real
    distances::Vector{Real}
    abs_errors::Vector{Real}
    rel_errors::Vector{Real}
    vol_stds::Vector{Real}
    n_drills::Int
    r_massive::Real
    last_action::Symbol
end


@with_kw struct BetaZeroResults <: Results
    config::MEConfiguration
    B::Array # set of beliefs
    Π::Array # set of policy estimates
    Z::Array # set of returns
end


function Base.convert(::Type{Dict}, structure::Union{Results, Configuration, JobParameters})
    return Dict(fn=>begin
                        field = getfield(structure, fn)
                        structtype = [Results, Configuration, JobParameters]
                        if supertype(typeof(field)) in structtype
                            convert(Dict, field) # recurse on custom structures
                        else
                            field
                        end
                    end for fn in fieldnames(typeof(structure)))
end


# TODO: MFMS.configurations
function MixedFidelityModelSelection.configurations(::Type{MEConfiguration};
                        pomcpow_iterations=[10,100,1000],
                        grid_dims_xys=[10,30,50],
                        num_seeds=10,
                        mainbody_types=[BlobNode, EllipseNode, CircleNode],
                        params::MEJobParameters)    
    configs = MEConfiguration[]
    for pomcpow_iters in pomcpow_iterations
        for grid_dims_xy in grid_dims_xys
            for mainbody_type in mainbody_types
                for seed in 1:num_seeds
                    grid_dims = (grid_dims_xy, grid_dims_xy, 1)
                    config = MEConfiguration(seed, grid_dims, pomcpow_iters, mainbody_type, false, params)
                    push!(configs, config)
                end
            end
        end
    end
    return configs::Vector{<:Configuration}
end


# TODO: MFMS.initialize
function MixedFidelityModelSelection.initialize(config::Configuration)
    # Configuration specific option
    seed = config.seed
    grid_dims = config.grid_dims
    mainbody = config.mainbody_type(grid_dims=grid_dims)
    pomcpow_iters = config.pomcpow_iters
    use_mcts = config.use_mcts

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

    up = MEBeliefUpdater(m, 1000, 2.0)
    b0 = POMDPs.initialize_belief(up, ds0)

    exploration_coefficient = 100.0
    k_action = 2.0
    alpha_action = 0.25

    if use_mcts
        next_action = NextActionSampler()
        f_next_action(bmdp::BeliefMDP, b::MEBelief, h) = POMCPOW.next_action(next_action, bmdp.pomdp, b, h)
        belief_reward(pomdp::POMDP, b, a, bp) = mean(reward(pomdp, s, a) for s in particles(b))
        bmdp = BeliefMDP(m, up, belief_reward)
        estimate_value = (bmdp, b, d) -> 0.0 # mean(MineralExploration.extraction_reward(bmdp.pomdp, s) for s in particles(b))
        solver = DPWSolver(n_iterations=pomcpow_iters,
                           check_repeat_action=true,
                           exploration_constant=exploration_coefficient,
                           next_action=f_next_action,
                           k_action=k_action,
                           alpha_action=alpha_action,
                           tree_in_info=true,
                           estimate_value=estimate_value)
        planner = solve(solver, bmdp)
    elseif pomcpow_iters == -1
        planner = RandomPolicy(m; updater=up)
    else
        solver = POMCPOWSolver(tree_queries=pomcpow_iters,
                               check_repeat_obs=true,
                               check_repeat_act=true,
                               next_action=NextActionSampler(),
                               k_action,
                               alpha_action,
                               k_observation=2.0,
                               alpha_observation=0.1,
                               criterion=POMCPOW.MaxUCB(exploration_coefficient),
                               final_criterion=POMCPOW.MaxQ(),
                               # final_criterion=POMCPOW.MaxTries(),
                               estimate_value=0.0
                               # estimate_value=leaf_estimation
                               )
        planner = POMDPs.solve(solver, m)
    end
    trial = METrial(config, m, planner, s0, b0, up)
    return trial
end


# TODO: MFMS.evaluate
function MixedFidelityModelSelection.evaluate(trial::Trial; save_dir=nothing)
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
    timing = @timed begin
        trial_results = run_trial(m, up, planner, s0, b0, save_dir=save_dir, display_figs=false, verbose=false, collect_training_data=trial.config.use_mcts)
    end

    if trial.config.use_mcts
        data = trial_results[end]
        results = BetaZeroResults(config=trial.config,
                                  B=[d.b for d in data],
                                  Π=[d.π for d in data],
                                  Z=[d.z for d in data])
    else
        results = MEResults(config=trial.config,
                            timing=timing,
                            ore_map=s0.ore_map,
                            mass_map=s0.ore_map .>= m.massive_threshold,
                            massive_threshold=m.massive_threshold,
                            discounted_return=trial_results[1],
                            distances=trial_results[2],
                            abs_errors=trial_results[3],
                            rel_errors=trial_results[4],
                            vol_stds=trial_results[5],
                            n_drills=trial_results[6],
                            r_massive=trial_results[7],
                            last_action=trial_results[8])
    end
    return results::Results
end


function job(config::Configuration; results_dir=abspath("results"), save_dir=nothing)
    trial = initialize(config)
    results = MixedFidelityModelSelection.evaluate(trial; save_dir=save_dir)
    !isnothing(results_dir) && save(results, results_dir)
    return results::Results
end


mainbody_type_string(config::MEConfiguration) = replace(lowercase(string(config.mainbody_type.name.name)), "node"=>"")

fileformat(results::Results) = fileformat(results.config)
function fileformat(config::MEConfiguration)
    seed = string("seed", config.seed)
    mainbody_type = mainbody_type_string(config)
    grid_dims = string(config.grid_dims[1], "x", config.grid_dims[2])
    planner_str = config.use_mcts ? "mcts" : "pomcpow"
    pomcpow_iters = string(planner_str, config.pomcpow_iters)
    # job_params_hash = string("params", "0x", string(hash(config.params); base=16))

    filename = join(["results", seed, mainbody_type, grid_dims, pomcpow_iters, config.params.name], "_")
    return filename
end

results_filename(results_or_config::Union{Results,MEConfiguration}, results_dir::String) = joinpath(results_dir, string(fileformat(results_or_config), ".bson"))

# TODO: MFMS.save
function MixedFidelityModelSelection.save(results::Results, results_dir::String)
    !isdir(results_dir) && mkdir(results_dir)
    filename = results_filename(results, results_dir)
    res = convert(Dict, results)
    @save filename res
end
