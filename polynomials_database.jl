using DataFrames , Polynomials , CSV , Random , Combinatorics 
function de_casteljeu_alg(coeffs)
    l=length(coeffs)
    coeffs_copy=copy(coeffs)
    left=[coeffs_copy[1]]
    right=[coeffs_copy[l]]
    for i= 1:l-1
        for j=1:l-i
            coeffs_copy[j]=(coeffs_copy[j]+coeffs_copy[j+1])/2

        end
        push!(left,coeffs_copy[1])
        pushfirst!(right, coeffs_copy[l-i])
    end
    return left,  right 
end
function root_finder(coeffs,a,b)
    eps=1e-3
    if   all(coeffs.>0) || all(coeffs.<0) 
        return Float64[] 
    end
    if (b-a)<eps
        return  [(a+b)/2]
    end
    left, right=de_casteljeu_alg(coeffs)
    m=(a+b)/2
    return vcat(root_finder(left,a,m) ,root_finder(right,m,b))
end
function construct_ds(n, deg)
    df=DataFrame(Polynomial=Vector{Vector{Float64}}(), roots=Vector{Vector{Float64}}())
    polyn_r=0.#Πολυώνυμα με ρίζα
    polyn_nr=0#Πολυώνυμα χωρίς ρίζα
    while polyn_r<n/2 || polyn_nr<n/2
        coeffs=rand(deg+1)*2 .- 1
        rf = root_finder(coeffs, 0, 1)
        if !isempty(rf) && polyn_r>= n/2
            continue
        end
        if isempty(rf) && polyn_nr>=n/2
            continue
        end 
        push!(df, (Polynomial=coeffs, roots=rf))
        if !isempty(rf)
            polyn_r+=1  
        else
            polyn_nr+=1
        end
    end
    return df
end
df=construct_ds(10000, 5)
CSV.write("C:\\Users\\USER\\Desktop\\bernstein_data_base.csv", df)


