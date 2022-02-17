import BSON: @save

subdirectory = "../results/hf_truth_mvb/"
!isdir(subdirectory) && mkdir(subdirectory)

grid_dim_xys = [10, 30, 50]
pomcpow_iters = [10, 100, 1000]

@time begin
    for mainbody_type in [BlobNode, EllipseNode, CircleNode]
        true_mainbody_type = BlobNode
        @show (true_mainbody_type, mainbody_type)

        local sweeps = [(xy, iter) for xy in grid_dim_xys for iter in pomcpow_iters]
        local num_processes = length(sweeps)

        # DEBUG.
        # TODO: index error (NOTE seed = i*10)
        seeds_fn = i -> [i*10, i*100, i*1000, i*10000, i*100000]
        seed_str = "seeds$(length(seeds_fn(0)))"
        global poutput = @time pmap(i->begin
                                        combined_results = []
                                        for seed in seeds_fn(i)
                                            results = run_single_trial(seed, sweeps[i][1], sweeps[i][2]; TrueMainbodyType=true_mainbody_type, MainbodyType=mainbody_type)
                                            push!(combined_results, results)
                                        end
                                        [[combined_results[i][j] for i in 1:length(combined_results)] for j in 1:length(combined_results[1])]
                                    end, 1:num_processes)

        local mainbody_type_str = lowercase(string(mainbody_type)[1:end-4])
        local true_mainbody_type_str = lowercase(string(true_mainbody_type)[1:end-4])
        global filename = joinpath(subdirectory, string("analysis_", true_mainbody_type_str, "_", mainbody_type_str, "_", seed_str, "_sweep_xy.10.30.50_pomcp.10.100.1000.bson"))
        include("save_sweep.jl")
    end
end
