global RESULTS_DIR = abspath(joinpath(@__DIR__, "..", "results"))

function configurations_fixed_bores()
    name = "fixed_bores"
    params = MEJobParameters(name=name, min_bores=20, max_bores=20)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=10,
                             pomcpow_iterations=[100],
                             grid_dims_xys=[30,50])
end


function configurations_fixed_bores_500()
    name = "fixed_bores_500"
    params = MEJobParameters(name=name, min_bores=20, max_bores=20)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=500,
                             pomcpow_iterations=[1000],
                             grid_dims_xys=[30,50])
end


function configurations_regret()
    name = "regret"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=5,
                             pomcpow_iterations=[10,100,1000],
                             grid_dims_xys=[10,30,50])
end


function configurations_regret20()
    name = "regret20"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=20,
                             pomcpow_iterations=[10,100,1000],
                             grid_dims_xys=[10,30,50])
end

function configurations_regret100()
    name = "regret100"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=100,
                             pomcpow_iterations=[10,100,1000],
                             grid_dims_xys=[10,30,50])
end

function configurations_biasvar()
    name = "biasvar"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=100,
                             # num_realizations=20,
                             pomcpow_iterations=[10,100,1000],
                             grid_dims_xys=[10,30,50])
end

function configurations_10K_blobbiasfix()
    name = "10K_blobbiasfix"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=100,
                             pomcpow_iterations=[100,1000,10000],
                             grid_dims_xys=[10,30,50])
end

function configurations_500seeds()
    name = "500seeds"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=500,
                             pomcpow_iterations=[100,1000,10000],
                             grid_dims_xys=[10,30,50])
end


function configurations_random_policy()
    name = "random_policy"
    params = MEJobParameters(name=name, min_bores=20, max_bores=20)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=500,
                             pomcpow_iterations=[-1],
                             grid_dims_xys=[50])
end


function configurations_random_policy_100()
    name = "random_policy_100"
    params = MEJobParameters(name=name, min_bores=100, max_bores=100)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=100,
                             pomcpow_iterations=[-1],
                             grid_dims_xys=[50])
end


function configurations_blob_points_recomputed()
    name = "blob_points_recomputed"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=500,
                             pomcpow_iterations=[100,1000,10000],
                             grid_dims_xys=[10,30,50])
end


function configurations_blob_copied()
    name = "blob_copied"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=500,
                             pomcpow_iterations=[100,1000,10000],
                             grid_dims_xys=[10,30,50])
end


function configurations_blob_clamped()
    name = "blob_clamped"
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=500,
                             pomcpow_iterations=[100,1000,10000],
                             grid_dims_xys=[10,30,50])
end


function configurations_standard(name)
    params = MEJobParameters(name=name)
    configs = configurations(MEConfiguration;
                             params=params,
                             num_seeds=500,
                             pomcpow_iterations=[100,1000,10000],
                             grid_dims_xys=[10,30,50])
end


function makebatches(configs::Vector{<:Configuration}, n)
    batchsizes = round(Int, length(configs)/n)
    return collect(Iterators.partition(configs, batchsizes))
end


function kickoff(configuration_fn::Function, nbatches; results_dir=RESULTS_DIR, rerun=false)
    configs = configuration_fn()
    batches = makebatches(configs, nbatches)

    if nprocs() < nbatches
        addprocs(nbatches - nprocs())
    end
    @info "Number of processes: $(nprocs())"

    progress = Progress(length(configs))
    channel = RemoteChannel(()->Channel{Bool}(), 1)

    @async while take!(channel)
        next!(progress)
    end

    @time pmap(batch->begin
                for config in batch
                    # skip if already ran (potentially from a cut-off previous run)
                    if !isfile(results_filename(config, results_dir)) || rerun
                        job(config; results_dir=results_dir)
                    end
                    put!(channel, true) # trigger progress bar update
                end
            end, batches)
    put!(channel, false) # tell printing task to finish

    return reduce_results(configuration_fn; results_dir=results_dir)
end


function kickoff(configuration_fn::Function; results_dir=RESULTS_DIR, rerun=false)
    configs = configuration_fn()

    @time @showprogress for config in configs
        # skip if already ran (potentially from a cut-off previous run)
        if !isfile(results_filename(config, results_dir)) || rerun
            @info "———"
            job(config; results_dir=results_dir)
        end
    end

    return reduce_results(configuration_fn; results_dir=results_dir)
end


function reduce_results(configuration_fn::Function; results_dir=RESULTS_DIR, partial=false)
    configs = configuration_fn()
    results = Dict()
    use_mcts = false
    @showprogress for config in configs
        grid_dims = config.grid_dims
        mainbody_type = Symbol(mainbody_type_string(config))
        pomcpow_iters = config.pomcpow_iters
        key = (mainbody_type, grid_dims, pomcpow_iters)
        filename = results_filename(config, results_dir)
        use_mcts = config.use_mcts
        try
            res = BSON.load(filename, @__MODULE__)[:res]
            res[:config][:mainbody_type] = string(res[:config][:mainbody_type].name.name) # Remove dependence on MineralExploration
            res[:seed] = config.seed # Note, adding seed value to results.
            if !haskey(results, key)
                # Create empty dictionary with key=>[] for all keys
                # Handles appending to an array of all results (handles array of arrays)
                results[key] = Dict(zip(keys(res), [[] for _ in 1:length(keys(res))]))
            end
            merge!((a,b)->[a...,b], results[key], res)
        catch err
            if partial
                @warn err
            else
                rethrow()
            end
        end
    end
    name = first(configs).params.name

    if use_mcts
        X_train, Y_train = combine_betazero_training_data(results)
        data = Dict(:X=>X_train, :Y=>Y_train)
        betazero_results_filename = joinpath(results_dir, "betazero_training_$name.bson")
        @info "Saving results: $betazero_results_filename"
        @save betazero_results_filename data
    else
        combined_results_filename = joinpath(results_dir, "results_$name.bson")
        @info "Saving results: $combined_results_filename"
        @save combined_results_filename results
    end

    return results
end
