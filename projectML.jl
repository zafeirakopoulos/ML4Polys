using Random, Combinatorics, Flux, Polynomials, JSON
using Flux: Chain, Dense, ADAM, mse, train!, onehotbatch, logitcrossentropy
# This is the first version of the code for the training set 
#function train_set(n) # n is the number of polynomials
    #train_set = []
    #for i=1:n 
        #deg = floor(10*rand())
        #coef = []                       
        #for j=1:(deg+1)
            #a_j = floor(10*rand())
            #push!(coef,a_j)
        #end
        #push!(train_set,coef)
    #end    
    #return train_set
#end


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


function write_coeffs_list_in_polynomial_form(coeffs, language) # this function takes a list of coefficients and returns a 
                                                                # string representation of the polynomial in the form 
                                                                # a_n*x^n + ... + a_1*x + a_0           
    n = length(coeffs) 
    terms = String[]
    for i in 1:n
        if language == "julia"
            if i == 1 
                push!(terms, "$(coeffs[i])")
            elseif i ==2 
                push!(terms, "$(coeffs[i])*x")
            else
                push!(terms, "$(coeffs[i])*x^$(i-1)")    
            end
        end

        if language == "python"
             if i == 1 
                push!(terms, "$(coeffs[i])")
            elseif i ==2 
                push!(terms, "$(coeffs[i])*x")
            else
                push!(terms, "$(coeffs[i])*x**$(i-1)")    
            end
        end 
        
        if language == "Tex"
            if i == 1 
                push!(terms, "$(coeffs[i])")
            elseif i ==2 
                push!(terms, "$(coeffs[i])*x")
            else
                push!(terms, "$(coeffs[i])*x^$(i-1)")    
            end
        end
    end 
    return join(reverse(terms), "+")
end                 


function write_data_base_to_json(filename, data_set)
    data = []
    reps = 0
    M = change_of_basis(data_set[1]) 
    if filename == "data_base.json"
        for poly in data_set
            reps += 1
            bernstein = M*poly 
            r = roots(Polynomials.Polynomial(reverse(poly)))
            entry = Dict(
                "name" => "Polynomial number $reps",  
                "Tex" => write_coeffs_list_in_polynomial_form(poly, "Tex"),
                "julia" => write_coeffs_list_in_polynomial_form(poly, "julia"),
                "python" => write_coeffs_list_in_polynomial_form(poly, "python"),
                "monomial_basis" => poly,
                "bernstein_basis" => bernstein,
                "degree" => length(poly) - 1,
                "bitsize" => log2(abs(poly[length(poly)]) + 1),
                "#real_roots" => count(x -> isreal(x), r),
                "#integer_roots" => count(x -> isinteger(x), r),
                "roots" => [r[i] for i in 1:length(r)],
                "sign_changes_monomial_basis" => [number_of_sign_changes(poly)],
                "sign_changes_bernstein_basis" => [number_of_sign_changes(bernstein)]
            )
            push!(data, entry)
        end 
        open(filename, "w") do file
            JSON.print(file, data, 4)  
        end
    end     
    if filename == "D_degree_bitsize.json"
        for poly in data_set
            bernstein = M*poly 
            r = Sturm_theorem([poly])
            entry = Dict(
                "julia" => write_coeffs_list_in_polynomial_form(poly, "julia"),
                "monomial_basis" => poly,
                "bernstein_basis" => bernstein,
                "degree" => length(poly) - 1,
                "bitsize" => log2(abs(poly[length(poly)]) + 1),
                "#real_roots " => r[1]
            )
            push!(data, entry)
        end
        open(filename, "w") do file
            JSON.print(file, data, 4)  
        end
    end     
end 


function read_json_data_set() 
        data = JSON.parsefile(" D_degree_bitsize.json")
        monomial = [d["monomial_basis"] for d in data]
        bernstein = [d["bernstein_basis"] for d in data]
    return monomial, bernstein
end
    


function Cauchy_bound(data_set=read_json_data_set("check_set_monomial_basis.json")) # this function computes the Cauchy bound for the roots of the polynomials in the data set
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


function horner_scheme(p, x) # this function evaluates the polynomial p at x using Horner's Scheme
    y = BigFloat(0)
    for i in length(p):-1:1
        y = p[i] + x*y
    end
    return y
end

function evaluate_sign(p, x) # this function evaluates the sign of the polynomial p at x  
    sign_ = nothing 
    if horner_scheme(p, x) > 0 
        sign_ = 1
    elseif horner_scheme(p, x) < 0 
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





function Strum_theorem(data_set=read_json_data_set("check_set_monomial_basis.json"))   # this fundtion uses the Strum sequence to compute the number of real roots of the polynomials 
                                    # in the data set, using the Cauchy bound as the interval, to count 
                                    # the number of sign changes in the Strum sequence
    strum_seq = [] # this is the list which contains the Strum sequence for each polynomial in the data set
    for i in data_set 
        p0 = Polynomials.Polynomial(reverse(i)) # p0 is the original polynomial
        p1 = derivative(p0) # p1 is the derivative of p0
        seq = [p0,p1]  # the Strum sequence for the i-th polynomial in the data set  
        while true 
            q, p = divrem(seq[end-1], seq[end]) # p is the remainder of the division of the last two polynomials in the sequence
            p = -p
            if iszero(p)
                break 
            else 
                push!(seq,p)
            end    
        end
        push!(strum_seq,seq)
    end
    Bound = Cauchy_bound(data_set)
    number_of_real_roots_seq = [] # this list contains the number of real roots for each one of the polynomials
                                  # in the data set
                                  # The number of real roots of the i-th polynomial in the data set 
                                  # is the difference between the number of sign changes in the Strum sequence at a and b, 
                                  # where a and b are the endpoints of the interval defined by the Cauchy bound
    
    for i in 1:length(data_set) # choose the i-th polynomial in the data set  
        seq = strum_seq[i] # choose the Strum sequence of the i-th polynomial in the data set
        sign_a = []
        sign_b = []
        a = -Bound[i]
        b = Bound[i]
        for j in seq # choose the j-th polynomial in the Strum sequence of the i-th polynomial in the data set
            push!(sign_a, evaluate_sign(j,a))
            push!(sign_b, evaluate_sign(j,b))
        end 
        
        V_a = number_of_sign_changes(sign_a) # the number of sign changes at a
        V_b = number_of_sign_changes(sign_b) # the number of sign changes at b 
        
        push!(number_of_real_roots_seq, V_a - V_b) # the number of real roots is the difference between 
                                                    # the number of sign changes at a and b   

    end 
    return number_of_real_roots_seq
end                            


model = Flux.Chain(
    Flux.Dense(101, 50, relu),
    Flux.Dense(50, 101)
)


x_raw = read_json_data_set()
y_raw = Strum_theorem(x_raw)

x = Float32.(hcat(x_raw...)) ./ 1_000_000
y = onehotbatch(y_raw, 0:100)
data = [(x, y)]

loss(m,x,y) = Flux.logitcrossentropy(m(x), y)
opt = Flux.setup(Flux.Adam(1e-1), model)

for epoch in 1:500
    Flux.train!(loss, model, data, opt)
    println("epoch = $epoch, loss = ",loss(model, x, y))
end    
