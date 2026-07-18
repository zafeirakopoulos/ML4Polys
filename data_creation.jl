using Random, Combinatorics, Polynomials, JSON, GenericLinearAlgebra, DataStructures, AMRVW

function sparse_polynomial_generator_in_monomial_basis(deg::Number, n::Number) # n is the number of polynomials, deg is the degree of the polynomials
    list_of_polynomials = Vector{Vector{BigFloat}}() # this is the list which contains the data set of polynomials in monomial basis
    for i in 1:n
        coeff= Vector{BigFloat}()
        for j in 1:deg
            if rand(Bool) # rand(Bool) is a boolean variable which decides whether to include the monomial x^j
                push!(coeff, BigFloat(rand(-1000000:1000000)))
            else 
                push!(coeff, BigFloat(0))
            end 
        end
        push!(coeff, BigFloat(rand(-1000000:1000000))) # leading coefficient is non-zero
        push!(list_of_polynomials, coeff)
    end
    return list_of_polynomials
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

 function Cauchy_bound(polynomial::Vector{<:Number}) # this function computes the Cauchy bound for the roots of the given polynomial
    max_coeff = BigFloat(0)
    for i in 1:length(polynomial)-1
        max_coeff = BigFloat(max(max_coeff, abs(polynomial[i])))
    end
    leading_coeff = polynomial[end]
    R = BigFloat(1) + (BigFloat(max_coeff) / BigFloat(abs(leading_coeff)))
    return R
end

function evaluate_sign(p::Polynomials.Polynomial, x::Number) # this function evaluates the sign of the polynomial p at x  
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

function number_of_sign_changes(v::Vector{<:Number}) # this function computes the number of sign changes in a vector of signs (1, -1, 0) 
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

function Sturm_sequence(polynomial::Vector{<:Number}) # this function computes the Sturm sequence of a polynomial
    length(polynomial) == 1 && return [Polynomials.Polynomial(polynomial)] # if the polynomial is a constant, 
                                                                           # then the Sturm sequence is the polynomial itself 
    p0 = Polynomials.Polynomial(polynomial) # p0 is the original polynomial
    p1 = derivative(p0) # p1 is the derivative of p0
    p0_new = div(p0, gcd(p0, p1)) # this makes sure that the polynomial is square-free
    p1_new = derivative(p0_new) # p1 is the derivative of the square-free p0 
    seq = [p0_new, p1_new]  # the Sturm sequence for polynomial   
    while true 
        _, p = divrem(seq[end-1], seq[end]) # p is the remainder of the division of the last two polynomials in the sequence
        p = -p
        if iszero(p)
            break 
        else 
            push!(seq,p)
        end    
    end
    return seq 
end 

function calculate_bitsize(coeffs::Vector{<:Number})
    bitsize = maximum(iszero(c) ? BigInt(1) : floor(BigInt, log2(abs(c))) + 1 for c in coeffs)
    return bitsize
end     

function write_coeffs_list_in_polynomial_form(coeffs::Vector{<:Number}, language::String) 
                                                                # this function takes a list of coefficients and returns a 
                                                                # string representation of the polynomial in the form 
                                                                # a_n*x^n + ... + a_1*x + a_0           
    n = length(coeffs) 
    terms = String[]
    for i in 1:n
        iszero(coeffs[i]) && continue 

        if coeffs[i] == 1
            term = ""
        elseif coeffs[i] == -1 
            term = "-"
        else 
            term = "$(coeffs[i])" * "*"       
        end

        if language == "julia" || language == "Tex"
            if i == 1 
                push!(terms, "$(coeffs[i])")
            elseif i ==2 
                push!(terms, term * "x")
            else
                push!(terms, term * "x" * "^" * "$(i-1)")    
            end
        end

        if language == "python"
            if i == 1 
                push!(terms, "$(coeffs[i])")
            elseif i ==2 
                push!(terms, term * "x")
            else
                push!(terms, term * "x" * "**" * "$(i-1)")    
            end
        end 
        
    end 
    
    return join(reverse(terms), " ")
end                 
    
function create_data_base(filename::String, list_of_polynomials::Vector{<:Vector{<:Number}})
    data = []
    reps = 0
    
    for poly in list_of_polynomials
        reps += 1
        M = change_of_basis(poly) 
        bernstein = M*poly 
        roots = AMRVW.roots(poly)
        bound_mon = Cauchy_bound(poly)
        a_mon = -bound_mon
        b_mon = bound_mon
        bound_bern = Cauchy_bound(bernstein)
        a_bern = -bound_bern
        b_bern = bound_bern 
        seq_mon = Sturm_sequence(poly)
        seq_bern = Sturm_sequence(bernstein)
        sign_a_mon = [evaluate_sign(Polynomials.Polynomial(p), a_mon) for p in seq_mon]
        sign_b_mon = [evaluate_sign(Polynomials.Polynomial(p), b_mon) for p in seq_mon] 
        sign_a_bern = [evaluate_sign(Polynomials.Polynomial(p), a_bern) for p in seq_bern]
        sign_b_bern = [evaluate_sign(Polynomials.Polynomial(p), b_bern) for p in seq_bern]
        number_of_real_roots = count(x -> isreal(x), roots)
        number_of_integer_roots = count(x -> isinteger(x), roots)
        tex = write_coeffs_list_in_polynomial_form(poly, "Tex")
        julia = write_coeffs_list_in_polynomial_form(poly, "julia")
        python = write_coeffs_list_in_polynomial_form(poly, "python")
        entry = OrderedDict(
            "name" => "Polynomial number $reps",  
            "Tex" => tex,
            "julia" => julia,
            "python" => python,
            "monomial_basis" => poly,
            "bernstein_basis" => bernstein,
            "degree" => length(poly) - 1,
            "bitsize" => calculate_bitsize(poly),
            "#real_roots" => number_of_real_roots,
            "#integer_roots" => number_of_integer_roots,
            "roots" => [roots[i] for i in 1:length(roots)],
            "sign_changes_monomial_basis" => [(number_of_sign_changes(sign_a_mon), number_of_sign_changes(sign_b_mon))],
            "sign_changes_bernstein_basis" => [(number_of_sign_changes(sign_a_bern), number_of_sign_changes(sign_b_bern))]
        )
        push!(data, entry)
    end 
    open(filename, "w") do file
        JSON.print(file, data, 4)  
    end 
end 

function create_data_set_from_data_base(extraction_file::String, insertion_file::String, target_degree::Number, target_bitsize::Number, target_size::Number)
    data_base = JSON.parsefile(extraction_file)
    extracted_data = OrderedDict(
        "julia" => [],
        "monomial_basis" => [],
        "bernstein_basis" => [],
        "degree" => [],
        "bitsize" => [],
        "#real_roots" => []
    )
    
    reps = 0
    for poly in data_base
        if poly["degree"] == target_degree && poly["bitsize"] == target_bitsize
            reps +=1
            if reps > target_size
                break
            end     
            push!(extracted_data["julia"], poly["julia"])
            push!(extracted_data["monomial_basis"], poly["monomial_basis"])
            push!(extracted_data["bernstein_basis"], poly["bernstein_basis"])
            push!(extracted_data["degree"], poly["degree"])
            push!(extracted_data["bitsize"], poly["bitsize"])
            push!(extracted_data["#real_roots"], poly["#real_roots"])
        end
    end
    
    open(insertion_file, "w") do file
        JSON.print(file, extracted_data, 4)
    end
end

function add_to_data_base(filename::String, polynomials_to_be_added::Vector{<:Vector{<:Number}})
    data_base = JSON.parsefile(filename)
    new_entries = []
    reps = length(data_base) # this is the number of polynomials already in the data base
    for poly in polynomials_to_be_added
        reps += 1
        M = change_of_basis(poly)
        bernstein = M*poly 
        roots = AMRVW.roots(poly)
        bound_mon = Cauchy_bound(poly)
        a_mon = -bound_mon
        b_mon = bound_mon
        bound_bern = Cauchy_bound(bernstein)
        a_bern = -bound_bern
        b_bern = bound_bern 
        seq_mon = Sturm_sequence(poly)
        seq_bern = Sturm_sequence(bernstein)
        sign_a_mon = [evaluate_sign(Polynomials.Polynomial(p), a_mon) for p in seq_mon]
        sign_b_mon = [evaluate_sign(Polynomials.Polynomial(p), b_mon) for p in seq_mon] 
        sign_a_bern = [evaluate_sign(Polynomials.Polynomial(p), a_bern) for p in seq_bern]
        sign_b_bern = [evaluate_sign(Polynomials.Polynomial(p), b_bern) for p in seq_bern]
        number_of_real_roots = count(x -> isreal(x), roots)
        number_of_integer_roots = count(x -> isinteger(x), roots)
        tex = write_coeffs_list_in_polynomial_form(poly, "Tex")
        julia = write_coeffs_list_in_polynomial_form(poly, "julia")
        python = write_coeffs_list_in_polynomial_form(poly, "python")
        entry = OrderedDict(
            "name" => "Polynomial number $reps",  
            "Tex" => tex,
            "julia" => julia,
            "python" => python,
            "monomial_basis" => poly,
            "bernstein_basis" => bernstein,
            "degree" => length(poly) - 1,
            "bitsize" => calculate_bitsize(poly),
            "#real_roots" => number_of_real_roots,
            "#integer_roots" => number_of_integer_roots,
            "roots" => [roots[i] for i in 1:length(roots)],
            "sign_changes_monomial_basis" => [(number_of_sign_changes(sign_a_mon), number_of_sign_changes(sign_b_mon))],
            "sign_changes_bernstein_basis" => [(number_of_sign_changes(sign_a_bern), number_of_sign_changes(sign_b_bern))]
        )
        push!(new_entries, entry)
    end 

    open(filename, "a") do file
        for entry in new_entries
            println(file, JSON.json(entry))   
         end 
    end      
end 