using Flux
using Flux: Data.DataLoader
using Flux: onehotbatch, onecold, crossentropy
using Flux: @epochs
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

# CNN model
function LeNet5(; imgsize=(30, 30, 5), out_dim=1)
    out_conv_size = (imgsize[1] ÷ 4 - 3, imgsize[2] ÷ 4 - 3, 16)

    return Chain(
        Conv((5, 5), imgsize[end]=>6, relu),
        MaxPool((2, 2)),
        Conv((5, 5), 6=>16, relu),
        MaxPool((2, 2)),
        Flux.flatten,
        Dense(prod(out_conv_size), 120, relu),
        Dense(120, 84, relu),
        Dense(84, out_dim)
    )
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
tt_split = 0.8 # train/test split (80/20)
data = BSON.load(joinpath(@__DIR__, "..", "..", "POMDP-value-estimation", "betazero_training_betazero.bson"))[:data]
x_data, y_data = data[:X], data[:Y]

# Remove channels 3,4 (skewness and kurtosis)
USE_ONLY_MEAN_STD = false
if USE_ONLY_MEAN_STD
    x_data = x_data[:, :, [1,2,5], :]
    imgsize = (30,30,3)
else
    imgsize = (30,30,5)
end

n_data = size(x_data,4)
n_train = Int(n_data ÷ (1/tt_split))
n_valid = n_data - n_train

perm = randperm(n_data)
perm_train = perm[1:n_train]
perm_valid = perm[n_train+1:end]

x_train = x_data[:,:,:,perm_train]
y_train = y_data[:, perm_train]

x_valid = x_data[:,:,:,perm_valid]
y_valid = y_data[:, perm_valid]

# Put model/data onto GPU device
model = gpu(LeNet5(; imgsize))
x_train = gpu(x_train)
y_train = gpu(y_train)
x_valid = gpu(x_valid)
y_valid = gpu(y_valid)

# Create the full dataset (IMPORTANT: do after putting on GPU)
train_data = DataLoader((x_train, y_train), batchsize=1024, shuffle=true)

sqnorm(x) = sum(abs2, x)
penalty(λ=0) = λ*sum(sqnorm, Flux.params(model))
error(ŷ, y) = Flux.Losses.mse(ŷ, y)
loss(x, y) = Flux.Losses.mse(model(x), y) + penalty()

if false
# @time begin
    lr = 0.001 # learning rate
    opt = ADAM(lr)
    θ = Flux.params(model)
    update_frequency = 1

    number_epochs = 1000

    losses_train = []
    losses_valid = []
    @info "Beginning training $(size(x_train))"
    for e in 1:number_epochs
        for (x, y) in train_data
            _, back = Flux.pullback(() -> loss(x, y), θ)
            Flux.update!(opt, θ, back(1.0f0))
        end
        loss_train = loss(x_train, y_train)
        loss_valid = loss(x_valid, y_valid)
        push!(losses_train, loss_train)
        push!(losses_valid, loss_valid)
        if e % update_frequency == 0
            println("Epoch: ", e, " Loss Train: ", loss_train, " Loss Val: ", loss_valid)
            plot(xlims=(1, number_epochs), ylims=(0, 6000))
            plot!(1:e, losses_train, label="training")
            plot!(1:e, losses_valid, label="validation")
            display(plot!())
        end
    end
    
    # @epochs number_epochs Flux.train!(loss, θ, train_data, opt)
    # @show error(model(x_train), y_train)
    # @show error(model(x_valid), y_valid)
end