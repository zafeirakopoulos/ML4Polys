using Polynomials, GenericLinearAlgebra

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

function Sturm_theorem(list_of_polynomials::Vector{<:Vector{<:Number}}) # this fundtion uses the Sturm sequence 
                                                                        # to compute the number of real roots of the polynomials 
                                                                        # in list_of_polynomials, using the Cauchy bound as the interval, to count 
                                                                        # the number of sign changes in the Sturm sequence

    number_of_real_roots_seq = [] # this list contains the number of real roots for each one of the polynomials
                                  # in the list of polynomials.
                                  # The number of real roots of the i-th polynomial in that list 
                                  # is the difference between the number of sign changes in the Sturm sequence at a and b, 
                                  # where a and b are the endpoints of the interval defined by the Cauchy bound
   
    for poly in list_of_polynomials 
        seq = Sturm_sequence(poly) # compute the Sturm sequence for the i-th polynomial in the list of polynomials 
        
        R = Cauchy_bound(poly)
        a = -R 
        b = R 
        sign_a = [evaluate_sign(p, a) for p in seq]
        sign_b = [evaluate_sign(p, b) for p in seq]

        V_a = number_of_sign_changes(sign_a) # the number of sign changes at a
        V_b = number_of_sign_changes(sign_b) # the number of sign changes at b 
        
        push!(number_of_real_roots_seq, V_a - V_b) # the number of real roots is the difference between 
                                                   # the number of sign changes at a and b   
    end
    return number_of_real_roots_seq
end                            
