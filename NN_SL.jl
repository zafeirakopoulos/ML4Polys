using Random, Flux, JSON, GenericLinearAlgebra, DataStructures
using Flux: Chain, Dense, ADAM, mse, train!, onehotbatch, logitcrossentropy

Random.seed!(1234)

function read_json_data_set(filename::String) 
    data = JSON.parsefile(filename)
    monomial = [BigFloat.(v) for v in data["monomial_basis"]]  
    bernstein = [BigFloat.(v) for v in data["bernstein_basis"]]
    n_real_roots = data["#real_roots"]
    return monomial, bernstein, n_real_roots
end

# classification model 
println("Classification model")
x_raw_mon_class = read_json_data_set("D_100_20.json")[1]
y_raw_mon_class = read_json_data_set("D_100_20.json")[3]

x_mon_class = hcat(x_raw_mon_class...)
y_mon = onehotbatch(y_raw_mon_class, 0:100)

n = size(x_mon_class, 2)
idx = shuffle(1:n)
train_size = Int(0.8*n)
test_idx  = idx[train_size+1:end]

x_mon_train_class = x_mon_class[:, idx[1:train_size]]
x_mon_test_class = x_mon_class[:, idx[train_size+1:end]]

y_mon_train = y_mon[:, idx[1:train_size]]
y_mon_test  = y_mon[:, idx[train_size+1:end]]

train_data_class = [(x_mon_train_class, y_mon_train)]
test_data_class = [(x_mon_test_class, y_mon_test)]

model_class = Flux.Chain(
    Flux.Dense(size(x_mon_class, 1), 50, relu),
    Flux.Dense(50, 10, relu),
    Flux.Dense(10, 101)
)

loss_class(m,x_mon_train,y_mon_train) = Flux.logitcrossentropy(m(x_mon_train), y_mon_train)
opt_class = Flux.setup(Flux.Adam(1e-3), model_class)

for epoch in 1:500
    Flux.train!(loss_class, model_class, train_data_class, opt_class)
    println("epoch = $epoch, loss = ",loss_class(model_class, x_mon_train_class, y_mon_train))
end    

loss_class(m,x_mon_test,y_mon_test) = Flux.logitcrossentropy(m(x_mon_test), y_mon_test)
println("loss = ",loss_class(model_class, x_mon_test_class, y_mon_test))


#regression model 
println("Regression model")
x_raw_mon_reg = read_json_data_set("D_100_20.json")[1]
y_raw_mon_reg = read_json_data_set("D_100_20.json")[3] 

x_mon_reg = hcat(x_raw_mon_reg...)
y_mon_reg = reshape(y_raw_mon_reg, 1, :)

n = size(x_mon_reg, 2)
idx_reg = shuffle(1:n)
train_size_reg = Int(0.8*n)

x_mon_train_reg = x_mon_reg[:, idx_reg[1:train_size_reg]]
x_mon_test_reg = x_mon_reg[:, idx_reg[train_size_reg+1:end]]

y_mon_train_reg = y_mon_reg[:, idx_reg[1:train_size_reg]]
y_mon_test_reg = y_mon_reg[:, idx_reg[train_size_reg+1:end]]

train_data_reg = [(x_mon_train_reg, y_mon_train_reg)]
test_data_reg = [(x_mon_test_reg, y_mon_test_reg)]

model_reg = Flux.Chain(
    Flux.Dense(size(x_mon_reg, 1), 50, relu),
    Flux.Dense(50, 10, relu),
    Flux.Dense(10, 1)
)

loss_reg(m, x_mon_train_reg, y_mon_train_reg) = Flux.mse(m(x_mon_train_reg), y_mon_train_reg)
opt_reg = Flux.setup(Flux.Adam(1e-1), model_reg)


for epoch in 1:500
    Flux.train!(loss_reg, model_reg, train_data_reg, opt_reg)
    println("epoch = $epoch, loss = ",loss_reg(model_reg, x_mon_train_reg, y_mon_train_reg))
end

loss_reg(m, x_mon_test_reg, y_mon_test_reg) = Flux.mse(m(x_mon_test_reg), y_mon_test_reg)
println("loss = ",loss_reg(model_reg, x_mon_test_reg, y_mon_test_reg))