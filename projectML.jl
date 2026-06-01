using Random, DataFrames, Flux, Lux, CSV, Polynomials, JSON
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


function data_set_monomial_base(; deg=100, n=20) # n is the number of polynomials, deg is the degree of the polynomials
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

global data_set = data_set_monomial_base(deg=100, n=20) # this is the global variable which contains the 
                                                        # data set of polynomials in monomial basis 


function  save_data_set(data_set, filename::String)
    open(filename, "w") do io
        JSON.print(io, data_set)
    end
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

function Strum_sequence(data_set)   # this fundtion uses the Strum sequence to compute the number of real roots of the polynomials 
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
            if evalpoly(a,coeffs(j)) != 0
                push!(sign_a, sign(evalpoly(a,coeffs(j))))  # evaluate the polynomials of that sequence at a 
            else 
                push!(sign_a, 0) # if the polynomial is zero at a, we push 0 to the list of signs at a    
            end
            if evalpoly(b,coeffs(j)) != 0
                push!(sign_b, sign(evalpoly(b,coeffs(j))))  # evaluate the polynomials of that sequence at b
            else
                push!(sign_b, 0) # if the polynomial is zero at b, we push 0 to the list of signs at b
            end
            

        end 
        

        sign_changes_a = 0
        for j in 1: length(sign_a)-1
            k = j+1
            if (sign_a[j] == 0 )
                continue
            elseif (sign_a[j] != 0 && sign_a[k] == 0)
                while k <= length(sign_a) && sign_a[k] == 0
                    k += 1
                end
                if (k <= length(sign_a)) && (sign_a[j] != sign_a[k])
                    sign_changes_a += 1
                end 
            else 
                if (sign_a[j] != sign_a[k])
                    sign_changes_a += 1
                end        
            end
        end
              
        
        sign_changes_b = 0
        for j in 1: length(sign_b)-1
            k = j+1
            if (sign_b[j] == 0) 
                continue
            elseif (sign_b[j] != 0 && sign_b[k] == 0)
                while (k <= length(sign_b) && sign_b[k] == 0)
                    k += 1
                end
                if (k <= length(sign_b)) && (sign_b[j] != sign_b[k])
                    sign_changes_b += 1
                end 
            else 
                if (sign_b[j] != sign_b[k])
                    sign_changes_b += 1
                end        
            end
        end  
        
        push!(number_of_real_roots_seq, sign_changes_a - sign_changes_b) # the number of real roots is the difference between 
                                                                         # the number of sign changes at a and b   
    end 
    return number_of_real_roots_seq
end                            