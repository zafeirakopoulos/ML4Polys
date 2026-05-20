using DataFrames , Polynomials, Random , Combinatorics 
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
function discriminant(coeffs)
    coeffs_der=zeros(l)
    Sylv=zeros(2l-3,2l-3)
    for i=1:l+1
        coeffs_der[i]=(i-1)*coeffs[i]
    end
    for i=1:l-2
        for j=1:l-1
            if i+j-1<2l-3
                Sylv[i,j+i-1]=coeffs[j]
            end
        end
    end
    for i=1:l-1
        for j=1:l-1
            if i+j-1<2l-3
                Sylv[l-1+i,i+j+1]=coeffs_der[j]
            end
        end 
    end
    res=det(Sylv)
    Discr=(-1)^((l-1)*(l-2)/2)*res/coeffs[1]
    return Discr    
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
function real_root_mult(coeffs)
    multi=0
    nonzero=[i for i in coeffs if i != 0]
    for i=1:l-1
        if nonzero[i]*nonzero[i+1]<0
            multi=multi+1
        end
    end
    return multi
end    
function construct_ds(n,deg)
    ds=DataFrame(Polynomial=Vector{Vector{Float64}}(),BernsteinPol=Vector{Vector{Float64}}(), degree=Float64, roots=Vector{Vector{Float64}}())
    polyn_r=0.#Πολυώνυμα με ρίζα
    polyn_nr=0#Πολυώνυμα χωρίς ρίζα
    while polyn_r<n/2 || polyn_nr<n/2
        coeffs=rand(deg+1)*2 .- 1
        bern_coeffs=change_of_basis(coeffs)
    end
end



