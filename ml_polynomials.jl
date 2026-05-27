using Flux ,  CSV, DataFrames 
df= CSV.read("C:\\Users\\USER\\Desktop\\polynomials_database.csv", DataFrame)
function parse_vec(s)
    if s== "Float64[]"
        return []
    end
    s=strip(s, ['[',']'])
    nums=split(s,",")
    return parse.(Float64,strip.(nums))
end 
Roots_parsed=[parse_vec(s) for s in df.roots]
max_l=maximum(length.(Roots_parsed))
function pad_vec(v)
    padded=copy(v)
    while length(padded)< max_l
        push!(padded, NaN)
    end 
    return padded
end
Polynom = [parse_vec(s) for s in df.Polynomial]
Roots=[pad_vec(v) for v in Roots_parsed]
model=Chain(Dense(6,90,relu),Dense(90,45,relu),Dense(45,6))
opt = Flux.setup(Flux.Adam(0.01), model)
for epoch in 1: 100
    loss, grads= Flux.withgradient(model) do m
        Roots_hat=m(Polynom)
        Flux.mse(Roots_hat, Roots)
    end
Flux.update!(opt_state, model, grads[1])

end 


