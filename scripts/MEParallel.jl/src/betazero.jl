function sample_training_configurations(n::Int; grid_dims=(30,30,1), planning_iters=10, mainbody_type=BlobNode)
    params = MEJobParameters(name="betazero")
    configs = MEConfiguration[]
    for i in 1:n
        Random.seed!(i)
        config = MEConfiguration(i, grid_dims, planning_iters, mainbody_type, true, params)
        push!(configs, config)
    end
    return configs::Vector{<:Configuration}
end


function combine_betazero_training_data(results::Dict)
    k = collect(keys(results))
    @assert length(k) == 1 # for BetaZero run, we only use one key
    res = results[k[1]]
    B = res[:B]
    Z = res[:Z]
    b0 = B[1][1]
    N = sum(map(length, B))
    X = Array{Float32,4}(undef, size(b0)..., N)
    Y = Array{Float32,2}(undef, 1, N)
    i = 1
    for (bs, zs) in zip(B, Z)
        for (b,z) in zip(bs, zs)
            X[:,:,:,i] = b
            Y[1, i] = z
            i += 1
        end
    end
    return X, Y
end
