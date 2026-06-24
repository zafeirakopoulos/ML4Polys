using Random, Combinatorics, Flux, Polynomials, JSON, GenericLinearAlgebra, DataStructures, AMRVW
using Flux: Chain, Dense, ADAM, mse, train!, onehotbatch, logitcrossentropy

function data_set_monomial_basis(; deg=100, n=1000) # n is the number of polynomials, deg is the degree of the polynomials
    data_set = Vector{Vector{BigFloat}}() # this is the list which contains the data set of polynomials in monomial basis
    for i in 1:n
        coeff=BigFloat[]
        for j in 1:deg
            if rand(Bool) # rand(Bool) is a boolean variable which decides whether to include the monomial x^j
                push!(coeff, BigFloat(rand(-1000000:1000000)))
            else 
                push!(coeff, BigFloat(0))
            end 
        end
        push!(coeff, BigFloat(rand(-1000000:1000000))) # leading coefficient is non-zero
        push!(data_set, coeff)
    end
    return data_set
end


function change_of_basis(coeffs) # courtesy of Evangelia Symeonidi
                                 # added BigFloat to match my data type   
      l=length(coeffs)
      bern_basis=zeros(BigFloat, l, l)
      for i=1:l
          for j=i:l
              bern_basis[i , j]=BigFloat(binomial(big(l-1), j-1)*binomial(big(j-1), i-1)*(-1)^(j-i))
          end
      end 
     change_matrix= inv(bern_basis)
     return change_matrix
 end
function write_coeffs_list_in_polynomial_form(coeffs, language)
    n = length(coeffs)
    terms = String[]
    pow_sym = (language == "python") ? "**" : "^"
    for i in 1:n
        c = coeffs[i]
        if i == 1
            term = "$(abs(c))"
        elseif i == 2
            term = "$(abs(c))*x"
        else
            term = "$(abs(c))*x$(pow_sym)$(i-1)"
        end
        if c >= 0
            push!(terms, "+" * term)
        else
            push!(terms, "-" * term)
        end
    end
    final_str = join(reverse(terms))
    if startswith(final_str, "+")
        final_str = final_str[2:end]
    end
    return final_str
end

function write_data_base_to_json(filename, data_set)
    data = []
    reps = 0
    M = change_of_basis(data_set[1]) 
    if filename == "data_base.json"
        for poly in data_set
            reps += 1
            bernstein = M*poly 
            roots = AMRVW.roots(poly)
            entry = OrderedDict(
                "name" => "Polynomial number $reps",  
                "Tex" => write_coeffs_list_in_polynomial_form(poly, "Tex"),
                "julia" => write_coeffs_list_in_polynomial_form(poly, "julia"),
                "python" => write_coeffs_list_in_polynomial_form(poly, "python"),
                "monomial_basis" => poly,
                "bernstein_basis" => bernstein,
                "degree" => length(poly) - 1,
                "bitsize" => log2(abs(poly[length(poly)]) + 1),
                "#real_roots" => count(x -> isreal(x), roots),
                "#integer_roots" => count(x -> isinteger(x), roots),
                "roots" => [roots[i] for i in 1:length(roots)],
                "sign_changes_monomial_basis" => [number_of_sign_changes(poly)],
                "sign_changes_bernstein_basis" => [number_of_sign_changes(bernstein)]
            )
            push!(data, entry)
        end 
        open(filename, "w") do file
            JSON.print(file, data, 4)  
        end
    end     
    if filename == "data_set.json"
        for poly in data_set
            bernstein = M*poly 
            roots = Sturm_theorem([poly])
            entry = OrderedDict(
                "julia" => write_coeffs_list_in_polynomial_form(poly, "julia"),
                "monomial_basis" => poly,
                "bernstein_basis" => bernstein,
                "degree" => length(poly) - 1,
                "bitsize" => log2(abs(poly[length(poly)]) + 1),
                "#real_roots " => roots[1]
            )
            push!(data, entry)
        end
        open(filename, "w") do file
            JSON.print(file, data, 4)  
        end
    end     
end 


function read_json_data_set() 
        data = JSON.parsefile("data_set.json")
        monomial = [d["monomial_basis"] for d in data]
        bernstein = [d["bernstein_basis"] for d in data]
    return monomial, bernstein
end


function Cauchy_bound(data_set) # this function computes the Cauchy bound for the roots of the polynomials in the data set
    bound = BigFloat[]
    for i in 1:length(data_set)
        max_coeff = BigFloat(0)
        for j in 1:length(data_set[i])-1
            max_coeff = BigFloat(max(max_coeff, abs(data_set[i][j])))
        end
        leading_coeff = data_set[i][length(data_set[i])]
        R = BigFloat(1) + (BigFloat(max_coeff) / BigFloat(abs(leading_coeff)))
        push!(bound, R)
    end
    return bound
end


function evaluate_sign(p, x) # this function evaluates the sign of the polynomial p at x  
    sign_ = nothing 
    if p(x) > 0 
        sign_ = 1
    elseif p(x) < 0 
        sign_ = -1
    else 
        sign_ = 0
    end 
    return sign_ 
end                   


function number_of_sign_changes(v) # this function computes the number of sign changes in a vector of signs (1, -1, 0) 
                                        # while ignoring the zeros, since they do not contribute to the number of sign changes 
    sign_changes = 0 
    prev = nothing  
    i = 1
    while i <= length(v)
        if v[i] == 0
            i += 1
            continue 
        end 
        
        if prev !== nothing && v[i] != prev 
            sign_changes += 1
        end
        
        prev = v[i]
        i += 1

    end 
    return sign_changes
end         



function Sturm_theorem(data_set) # this fundtion uses the Sturm sequence to compute the number of real roots of the polynomials 
                                 # in the data set, using the Cauchy bound as the interval, to count 
                                 # the number of sign changes in the Sturm sequence
    sturm_seq = [] # this is the list which contains the Sturm sequence for each polynomial in the data set
    for i in data_set 
        p0 = Polynomials.Polynomial(i) # p0 is the original polynomial
        p1 = derivative(p0) # p1 is the derivative of p0
        p0_new = div(p0, gcd(p0, p1)) # this makes sure that the polynomial is square-free
        p1_new = derivative(p0_new) # p1 is the derivative of the square-free p0 
        seq = [p0_new, p1_new]  # the Sturm sequence for the i-th polynomial in the data set  
        while true 
            _, p = divrem(seq[end-1], seq[end]) # p is the remainder of the division of the last two polynomials in the sequence
            p = -p
            if iszero(p)
                break 
            else 
                push!(seq,p)
            end    
        end
        push!(sturm_seq,seq)
    end
    Bound = Cauchy_bound(data_set)
    number_of_real_roots_seq = [] # this list contains the number of real roots for each one of the polynomials
                                  # in the data set
                                  # The number of real roots of the i-th polynomial in the data set 
                                  # is the difference between the number of sign changes in the Sturm sequence at a and b, 
                                  # where a and b are the endpoints of the interval defined by the Cauchy bound
    
    for i in 1:length(data_set) # choose the i-th polynomial in the data set  
        seq = sturm_seq[i] # choose the Sturm sequence of the i-th polynomial in the data set
        a = -Bound[i]
        b = Bound[i]
        sign_a = [evaluate_sign(p, a) for p in seq]
        sign_b = [evaluate_sign(p, b) for p in seq]
        
        V_a = number_of_sign_changes(sign_a) # the number of sign changes at a
        V_b = number_of_sign_changes(sign_b) # the number of sign changes at b 
        
        push!(number_of_real_roots_seq, V_a - V_b) # the number of real roots is the difference between 
                                                    # the number of sign changes at a and b   

    end 
    return number_of_real_roots_seq
end                            


# classification model 
x_raw_mon_class = read_json_data_set()[1]
y_raw_mon_class = Sturm_theorem(x_raw_mon_class)

x_mon_class = hcat(x_raw_mon_class...)
y_mon = onehotbatch(y_raw_mon_class, 0:100)

n = size(x_mon_class, 2)
idx = shuffle(1:n)
train_size = Int(0.8*n)
test_idx  = idx[train_size+1:end]

x_mon_train = x_mon_class[:, idx[1:train_size]]
x_mon_test = x_mon_class[:, idx[train_size+1:end]]

y_train = y_mon[:, idx[1:train_size]]
y_test  = y_mon[:, idx[train_size+1:end]]

train_data_class = [(x_mon_train, y_train)]
test_data_class = [(x_mon_test, y_test)]

model_class = Flux.Chain(
    Flux.Dense(size(x_mon_class, 1), 10, relu),
    Flux.Dense(10, 101)
)

loss_class(m,x_mon_train,y_mon_train) = Flux.logitcrossentropy(m(x_mon_train), y_mon_train)
opt_class = Flux.setup(Flux.Adam(1e-1), model_class)

for epoch in 1:500
    Flux.train!(loss_class, model_class, train_data_class, opt_class)
    println("epoch = $epoch, loss = ",loss_class(model_class, x_mon_train, y_mon_train))
end    

loss_class(m,x_mon_test,y_mon_test) = Flux.logitcrossentropy(m(x_mon_test), y_mon_test)
opt_class = Flux.setup(Flux.Adam(1e-1), model_class)
Flux.train!(loss_class, model_class, test_data_class, opt_class)
println("epoch = $epoch, loss = ",loss_class(model_class, x_mon_test, y_mon_test))


#regression model 
x_raw_mon_reg = read_json_data_set()[1]
y_raw_mon_reg = Sturm_theorem(x_raw_mon_reg)

x_mon_reg = hcat(x_raw_mon_reg...)
y_mon_reg = reshape(y_raw_mon_reg, 1, :)

n = size(x_mon_reg, 2)
idx_reg = shuffle(1:n)
tain_size_reg = Int(0.8*n)

x_mon_train_reg, y_mon_train_reg = x_mon_reg[:, idx_reg[1:train_size_reg]], y_mon_reg[:, idx_reg[1:train_size_reg]]
x_mon_test_reg, y_mon_test_reg = x_mon_reg[:, idx_reg[train_size_reg+1:end]], y_mon_reg[:, idx_reg[train_size_reg+1:end]]

train_data_reg = [(x_mon_train_reg, y_mon_train_reg)]
test_data_reg = [(x_mon_test_reg, y_mon_test_reg)]

model_reg = Flux.Chain(
    Flux.Dense(size(x_mon_reg, 1), 10, relu),
    Flux.Dense(10, 1)
)

loss_reg(m, x_mon_train_reg, y_mon_train_reg) = Flux.mse(m(x_mon_train_reg), y_mon_train_reg)
opt_reg = Flux.setup(Flux.Adam(1e-1), model_reg)


for epoch in 1:500
    Flux.train!(loss_reg, model_reg, train_data_reg, opt_reg)
    println("epoch = $epoch, loss = ",loss_reg(model_reg, x_mon_train_reg, y_mon_train_reg))
end

loss_reg(m, x_mon_test_reg, y_mon_test_reg) = Flux.mse(m(x_mon_test_reg), y_mon_test_reg)
opt_reg = Flux.setup(Flux.Adam(1e-1), model_reg)
Flux.train!(loss_reg, model_reg, test_data_reg, opt_reg)
println("epoch = $epoch, loss = ",loss_reg(model_reg, x_mon_test_reg, y_mon_test_reg))
