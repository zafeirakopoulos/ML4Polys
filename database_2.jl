using JSON ,AbstractAlgebra , Random , Combinatorics 
l=length(coeffs)
function change_of_basis(coeffs)
    bern_basis=zeros(l, l)
    for i=1:l
        for j=i:l
            bern_basis[i , j]=binomial(l-1, j-1)*binomial(j-1, i-1)*(-1)^(j-i)
        end
    end 
    change_matrix= inv(bern_basis)
    return change_matrix
end
function mmbound(coeffs)
    deg=length(coeffs)-1
    norm=0
    for i in coeffs
        norm=norm+i^2
    end
    norm=sqrt(norm)
    b_mm=sqrt(3abs(discriminant(coeffs)))/deg^(deg/2+1)*norm^(deg-2)
    return b_mm
end

function root_finder(coeffs,a=nothing,b=nothing)
    if a === nothing || b === nothing
        c_bound=1+maximum(abs(coeffs[i]) for i=2:l)/abs(coeffs[1])
        a=-c_bound
        b=c_bound
    end
    if (b-a)<mmbound 
        return [(a+b)/2]
    end
    m=(a+b)/2
     roots=Float64[]
    p_a=sum(coeff[i]*a^(l-i) for i=1:l)
    p_m=sum(coeff[i]*m^(l-i) for i=1:l)
    p_b=sum(coeff[i]*b^(l-i) for i=1:l)
    if p_a*p_m<0
        append!(roots, root_finder(coeffs,a,m))
    end
    if p_m*p_b<0
        append!(roots, root_finder(coeffs,m,b))
    end
    return unique(roots)         
end
function sturm_seq(coeff_s)
    l=length(coeff_s)
    p_0=Polynomial(Rational{Int64}.(coeff_s))
    p_1=derivative(p_0)
    sturm_s=[p_0,p_1]
    while degree(sturm_s[end])>0
        r=rem(sturm_s[end-1],sturm_s[end])
        # r_clean = truncate(r, atol=1e-10)
        if r == 0
            break
        end
        
        # p_2 = -r_clean
         p_2=-r
        push!(sturm_s, p_2)
    end
    return sturm_s
end
function sign_changes(coeff_s, x)
    st=sturm_seq(coeff_s)
    signs=Int[]
    for p in st 
        p_float = Polynomial(Float64.(coeffs(p)))
        p_x = p_float(x)
        if p_x>0
            push!(signs,1)
        elseif p_x<0
            push!(signs,-1)
        end 
    end
    changes=0
    for i=1:(length(signs)-1)
        if signs[i] != signs[i+1] 
            changes+=1
        end
    end
    return changes

end
function discr_real_roots(coeff_s)
    l = length(coeff_s)
    v_b=sign_changes(coeff_s, Inf)
    v_s=sign_changes(coeff_s, -Inf)
    return v_s-v_b
end
function real_root_count(coeff_s)
    
end
coeffs_example = [0, 4, 0, -5, 0, 1]
counts=real_root_count(coeffs_example)
print(counts) 
function construct_json_db(n,deg)
     l=length(coeffs)
     database = Dict{String, Any}[]
     polyn_r=0.#Πολυώνυμα με ρίζα
     polyn_nr=0#Πολυώνυμα χωρίς ρίζα
     while polyn_r<n/2 || polyn_nr<n/2
          coeffs=rand(deg+1)*2 .- 1
          bern_coeffs=change_of_basis(coeffs)*coeffs
          root_multi=real_root_count(coeffs)
          if root_multi==0 
             polyn_nr=polyn_nr+1
             poly_entry = Dict(
                 "Polynomial"   => coeffs,
                 "BernsteinPol" => bern_coeffs,
                 "degree"       => deg,
                 "root_multi"   => root_multi,
                 "roots"        => Float64[]
             )
             push!(database, poly_entry)
          else
             polyn_r=polyn_r+1 
             roots=root_finder(coeffs)
             poly_entry = Dict(
                 "Polynomial"   => coeffs,
                 "BernsteinPol" => bern_coeffs,
                 "degree"       => deg,
                 "root_multi"   => root_multi,
                 "roots"        => roots
             )
             push!(database, poly_entry)
          end


     end
     return database
 end

