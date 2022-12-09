using Flux
using Flux: Data.DataLoader
using Flux: onehotbatch, onecold, crossentropy
using Flux: @epochs
using CUDA
using Statistics
using MLDatasets
using BSON
using Random
using Plots; default(fontfamily="Computer Modern", framestyle=:box)

function plot_image(X, i)
    img = X[:,:,1,i]
    return heatmap(rotl90(img), ratio=1, xlims=ylims(), c=:cividis)
end

Random.seed!(0) # determinism
CUDA.seed!(0)

USE_ONLY_MEAN_STD = false
NORMALIZE_TO_SIGN = true

# CNN model
function LeNet5(; imgsize=(30, 30, 5), out_dim=1, num_filters=6)
    out_conv_size = (imgsize[1] ÷ 4 - 3, imgsize[2] ÷ 4 - 3, 16)

    m = Chain(
        Conv((5, 5), imgsize[end]=>num_filters, relu),
        # MeanPool((2, 2)), # MaxPool((2, 2)),
        Conv((5, 5), num_filters=>16, relu),
        # MeanPool((2, 2)), # MaxPool((2, 2)),
        Flux.flatten,
        # Dense(prod(out_conv_size), 120, relu),
        Dense(7744, 120, relu),
        # Dropout(0.5),
        Dense(120, 84, relu),
        # Dropout(0.5),
        Dense(84, out_dim),
    )

    # if NORMALIZE_TO_SIGN
    #     m = Chain(m.layers..., Dense(out_dim, out_dim, tanh))
    # end

    return m
end


function AlphaZero(; imgsize=(30, 30, 5), out_dim=1, kernel_size=(3,3), num_filters=256, bnorm_momentum=0.6f0)
    out_conv_size = (imgsize[1] ÷ 4 - 3, imgsize[2] ÷ 4 - 3, 16)
    pad = kernel_size .÷ 2

    m = Chain(
        # Common    
        Conv(kernel_size, imgsize[end]=>num_filters, stride=1, pad=pad),
        BatchNorm(num_filters, relu, momentum=bnorm_momentum),

        # Value head
        Conv((1, 1), num_filters=>out_dim),
        BatchNorm(out_dim, relu, momentum=bnorm_momentum),
        Flux.flatten,
        Dense(imgsize[1]*imgsize[2]*out_dim, out_dim, relu),
    )

    if NORMALIZE_TO_SIGN
        m = Chain(m.layers..., Dense(out_dim, out_dim, tanh))
    end

    return m
end

    # model = Chain(
#     # 28x28 => 14x14
#     Conv((5, 5), 1=>8, pad=2, stride=2, relu),
#     # 14x14 => 7x7
#     Conv((3, 3), 8=>16, pad=1, stride=2, relu),
#     # 7x7 => 4x4
#     Conv((3, 3), 16=>32, pad=1, stride=2, relu),
#     # 4x4 => 2x2
#     Conv((3, 3), 32=>32, pad=1, stride=2, relu),

#     # Average pooling on each width x height feature map
#     GlobalMeanPool(),

#     # Fully-connected layer
#     Flux.flatten,
#     Dense(32, 10),

#     # Output layer
#     softmax
# )

# Load the data
data = BSON.load(joinpath(@__DIR__, "..", "..", "POMDP-value-estimation", "betazero_training_betazero_pomcpow.bson"))[:data]
x_data, y_data = data[:X], data[:Y]

# Remove channels 3,4 (skewness and kurtosis)
if USE_ONLY_MEAN_STD
    x_data = x_data[:, :, [1,2,5], :]
    imgsize = (30,30,3)
else
    imgsize = (30,30,5)
end

if NORMALIZE_TO_SIGN
    # y_data = sign.(y_data)
    mean_y = mean(y_data)
    std_y = std(y_data)
    y_data = (y_data .- mean_y) ./ std_y
end

n_data = 10_000 # size(x_data,4)
@show n_data
tt_split = 0.8 # train/test split (80/20)
n_train = Int(n_data ÷ (1/tt_split))
n_valid = n_data - n_train

perm = randperm(n_data)
perm_train = perm[1:n_train]
perm_valid = perm[n_train+1:n_data]

x_train = x_data[:,:,:,perm_train]
y_train = y_data[:, perm_train]

x_valid = x_data[:,:,:,perm_valid]
y_valid = y_data[:, perm_valid]

# Put model/data onto GPU device
x_train = gpu(x_train)
y_train = gpu(y_train)
x_valid = gpu(x_valid)
y_valid = gpu(y_valid)

# Create the full dataset (IMPORTANT: do after putting on GPU)
train_data = DataLoader((x_train, y_train), batchsize=512, shuffle=true)

RUN_TUNING = false
if RUN_TUNING
    tuning_results = Dict()
    LRs = [0.1, 0.01, 0.001, 0.005, 0.0001]
    λs = [0, 0.01, 0.005, 0.001, 0.0001]
    LOSSES = [true, false]
else
    tuning_results = BSON.load(joinpath(@__DIR__, "tuning_results2.bson"))[:tuning_results]
    LRs = [0.01]
    λs = [0.0001]
    LOSSES = [false]
end


results = Dict()

for lr in LRs
    for λ in λs
        for use_mse_loss in LOSSES
            global results
            @show lr, λ, use_mse_loss
            
            Random.seed!(0) # determinism
            CUDA.seed!(0)

            model = gpu(LeNet5(; imgsize))
            # model = gpu(AlphaZero(; imgsize))        

            sqnorm(x) = sum(abs2, x)
            penalty(λ=λ) = λ*sum(sqnorm, Flux.params(model))
            accuracy(x, y) = mean(sign.(model(x)) .== sign.(y))
            if use_mse_loss
                loss = (x, y)->Flux.Losses.mse(model(x), y) + penalty()
            else
                loss = (x, y)->Flux.Losses.mae(model(x), y) + penalty()
            end

            # if false
            @time begin
                opt = ADAM(lr)
                θ = Flux.params(model)
                update_frequency = 1
                plot_frequency = 10

                number_epochs = 1_000

                losses_train = []
                losses_valid = []
                accs_train = []
                accs_valid = []
                @info "Beginning training $(size(x_train))"
                for e in 1:number_epochs
                    for (x, y) in train_data
                        _, back = Flux.pullback(() -> loss(x, y), θ)
                        Flux.update!(opt, θ, back(1.0f0))
                    end
                    loss_train = loss(x_train, y_train)
                    loss_valid = loss(x_valid, y_valid)
                    acc_train = accuracy(x_train, y_train)
                    acc_valid = accuracy(x_valid, y_valid)
                    push!(losses_train, loss_train)
                    push!(losses_valid, loss_valid)
                    push!(accs_train, acc_train)
                    push!(accs_valid, acc_valid)
                    if e % update_frequency == 0
                        println("Epoch: ", e, " Loss Train: ", loss_train, " Loss Val: ", loss_valid, " | Acc. Train: ", acc_train, " Acc. Val: ", acc_valid)
                    end
                    if e % plot_frequency == 0
                        plot(xlims=(1, number_epochs), ylims=(0, NORMALIZE_TO_SIGN ? 1 : 2000))
                        plot!(1:e, losses_train, label="training")
                        plot!(1:e, losses_valid, label="validation")
                        display(plot!())
                    end
                end

                learning_curve = plot!()
                
                value_model = (cpu(model(x_valid))' .* std_y) .+ mean_y
                value_data = (cpu(y_valid)' .* std_y) .+ mean_y
                value_distribution = histogram(value_model, alpha=0.5, label="model", c=3)
                histogram!(value_data, alpha=0.5, label="data", c=4)
                display(plot!())

                if RUN_TUNING
                    tuning_results[(lr, λ, use_mse_loss)] = Dict(
                        "losses_train" => losses_train,
                        "losses_valid" => losses_valid,
                        "accs_train" => accs_train,
                        "accs_valid" => accs_valid,
                        "value_model" => value_model,
                        "value_data" => value_data,
                        "curve" => learning_curve,
                        "value_distribution" => value_distribution,                
                    )
                else
                    results = Dict(
                        "losses_train" => losses_train,
                        "losses_valid" => losses_valid,
                        "accs_train" => accs_train,
                        "accs_valid" => accs_valid,
                        "value_model" => value_model,
                        "value_data" => value_data,
                        "curve" => learning_curve,
                        "value_distribution" => value_distribution,                
                    )
                end
            end
        end
    end
end