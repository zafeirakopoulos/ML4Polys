using JSON, Random , Combinatorics 
using AbstractAlgebra
  function change_of_basis(coeffs)
      l=length(coeffs)
      bern_basis=zeros(l, l)
      for i=1:l
          for j=i:l
              bern_basis[i , j]=binomial(l-1, j-1)*binomial(j-1, i-1)*(-1)^(j-i)
          end
      end 
     change_matrix= inv(bern_basis)
     return change_matrix
 end
 R, x = polynomial_ring(QQ, "x")
function mmbound(coeff_s)
    p = R(big.(coeff_s))
    deg=big(length(coeff_s)-1)
    norm=big(0.0)
    for i in coeff_s
        norm=norm+big(i^2)
    end
    norm=sqrt(norm)
    b_mm=sqrt(3*Float64(abs((big(discriminant(p))))))/(deg^(deg/2+1)*norm^(deg-2))
    return b_mm
end

function root_finder(coeff_s,a=nothing,b=nothing)
    l=length(coeff_s)
    roots=Float64[]
    p = R(big.(coeff_s))
    p_prime = AbstractAlgebra.derivative(p)
    g = gcd(p, p_prime)
    q = divexact(p,g)  
    println("p = ", p)
    println("g = ", g)
    println("q = ", q)
    println("interval = [$a,$b]")
    if a === nothing || b === nothing
        c_bound=1+maximum(abs(coeff_s[i]) for i=1:l-1)/abs(coeff_s[end])
        a=Float64(-c_bound)
        b=Float64(c_bound)
    end
    q_a=Float64(evaluate(q,a))
    q_b=Float64(evaluate(q,b))
    if  (b-a)<mmbound(coeff_s) 
        if q_a*q_b<=0
            push!(roots, (a + b) / 2)
            return roots
        end
    end
    m=(a+b)/2
    q_a=Float64(evaluate(q,a))
    q_m=Float64(evaluate(q,m))
    q_b=Float64(evaluate(q,b))
    if iszero(q_m)
        return [m]
    end
    if q_a*q_m<=0 
        
        append!(roots, root_finder(coeff_s,a,m))
    end
    if q_m*q_b<=0 
        append!(roots, root_finder(coeff_s,m,b))
    end
    return unique(roots)         
end
#  function multi_roots(coeff_s)
#      roots=root_finder(coeff_s)
#      dv=[]
#     Θα φτιαξω μια λιστα η οποια θα περιεχει ολες τις παραγωγου της συναρτησης. Θα παιρνω καθε ριζα της συναρτησης
#     και θα ελεγχω μεχρι ποια παραγωγο μηδενιζει η ριζα. Ο αριθμος των παραγωγων που μηδενιζει ειναι ισος με την πολλαπλοτητα 
#     της ριζας 
#  end 
#  function descartes(coeffs)
#       l=length(coeffs)
#        pos_roots=0
#        nonzero=[i for i in coeffs if i != 0]
#        for i=1:l-1
#            if nonzero[i]*nonzero[i+1]<0
#                pos_roots=pos_roots+1
#            end
#        end
#        return pos_roots
#   end 
#  function sturm_seq(coeff_s)
#      l=length(coeff_s)
#      p_0=Polynomial(Rational{Int64}.(coeff_s))
#      p_1=derivative(p_0)
#      sturm_s=[p_0,p_1]
#      while degree(sturm_s[end])>0
#          r=rem(sturm_s[end-1],sturm_s[end])
#          # r_clean = truncate(r, atol=1e-10)
#          if r == 0
#              break
#          end
        
#          # p_2 = -r_clean
#           p_2=-r
#          push!(sturm_s, p_2)
#      end
#      return sturm_s
#  end
#  function sign_changes(coeff_s, x)
#      st=sturm_seq(coeff_s)
#      signs=Int[]
#      for p in st 
#          p_float = Polynomial(Float64.(coeffs(p)))
#          p_x = p_float(x)
#          if p_x>0
#              push!(signs,1)
#          elseif p_x<0
#              push!(signs,-1)
#          end 
#      end
#      changes=0
#      for i=1:(length(signs)-1)
#          if signs[i] != signs[i+1] 
#              changes+=1
#          end
#      end
#      return changes

#  end
#  function discr_real_roots(coeff_s)
#      l = length(coeff_s)
#      v_b=sign_changes(coeff_s, Inf)
#      v_s=sign_changes(coeff_s, -Inf)
#      return v_s-v_b
#  end
#  function real_root_count(coeff_s)
#      de_pos=descartes(coeff_s)
#      coeff_neg=[coeff_s[i]*(-1)^i for i=1:length(coeff_s)]
#      de_neg=descartes(coeff_neg)
#      if discr_real_roots(coeff_s)== de_pos+de_neg
#         return discr_real_roots(coeff_s)
#      else
#      end
#  end

coeffs_example = 
[ -1,  0,  8,  0,  27,  0,  48,  0,  42,  0, 0, 0,  -42,  0, -48, 0, -27,   0,  -8,   0,  1  ]
roo_t=root_finder(coeffs_example)

print(roo_t) 
println("mmbound = ", mmbound(coeffs_example))
#  function construct_json_db(n,deg)
#       l=length(coeffs)
#       database = Dict{String, Any}[]
#       polyn_r=0.#Πολυώνυμα με ρίζα
#       polyn_nr=0#Πολυώνυμα χωρίς ρίζα
#       while polyn_r<n/2 || polyn_nr<n/2
#            coeffs=rand(deg+1)*2 .- 1
#            bern_coeffs=change_of_basis(coeffs)*coeffs
#            root_multi=real_root_count(coeffs)
#            if root_multi==0 
#               polyn_nr=polyn_nr+1
#               poly_entry = Dict(
#                   "Polynomial"   => coeffs,
#                   "BernsteinPol" => bern_coeffs,
#                   "degree"       => deg,
#                   "root_multi"   => root_multi,
#                   "roots"        => Float64[]
#               )
#               push!(database, poly_entry)
#            else
#               polyn_r=polyn_r+1 
#               roots=root_finder(coeffs)
#               poly_entry = Dict(
#                   "Polynomial"   => coeffs,
#                   "BernsteinPol" => bern_coeffs,
#                   "degree"       => deg,
#                   "root_multi"   => root_multi,
#                  "roots"        => roots
#              )
#              push!(database, poly_entry)
#            end
#       end
#       return database
# end


