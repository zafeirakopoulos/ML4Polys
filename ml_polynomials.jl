using Flux ,  CSV, JSON 
using DataFrames
json_parsed= JSON.parsefile("C:\\Users\\USER\\Desktop\\polynomials_database.json")
df = DataFrame(json_parsed)
Roots_parsed=[Float32.(v) for v in df.roots]
Polynom=[Float32.(p) for p in df.Polynomial]
max_l = maximum(length.(Roots_parsed))
function pad_vec(v)
    padded=copy(v)
   while length(padded)< max_l
        push!(padded, 0.f0)
    end 
     return padded
end
Roots=[pad_vec(v) for v in Roots_parsed]
X_train=Float32.(reduce(hcat, Polynom))
Y_train=Float32.(reduce(hcat, Roots))
model=Chain(Dense(6,90,relu),Dense(90,45,relu),Dense(45,6))
opt = Flux.setup(Flux.Adam(0.01), model)
for epoch in 1: 100
    loss, grads= Flux.withgradient(model) do m
        Roots_hat=m(X_train)
        Flux.mse(Roots_hat, Y_train)
    end
    Flux.update!(opt, model, grads[1])
end 


