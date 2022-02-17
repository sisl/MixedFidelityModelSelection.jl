subdirectory = "../results/hf_truth_mvb/"
!isdir(subdirectory) && mkdir(subdirectory)

for pomcpow_iters in [100]
    # for grid_dims in [(10,10,1), (30,30,1), (50,50,1)]
    for grid_dims in [(30,30,1), (50,50,1)]
    # for grid_dims in [(10,10,1)]
        for mainbody_type in [BlobNode, EllipseNode, CircleNode]
            


            if grid_dims == (30,30,1) && mainbody_type == BlobNode
                @warn "ARE YOU SURE THIS SHOULD BE HERE?"
                continue
                # TODO: REMOVE
            end



            @show mainbody_type
            global poutput = @time pmap(i->run_single_trial(i, grid_dims, pomcpow_iters; MainbodyType=mainbody_type), 1:10)

            global filename = joinpath(subdirectory, string("analysis_", lowercase(string(mainbody_type)[1:end-4]), "_", grid_dims[1], "x", grid_dims[2], "_pomcp.", pomcpow_iters, ".bson"))
            include("save.jl")
        end
    end
end