using Dataframes, LinnearAlgebra, Polynomials, JSON

#We construct some functions that give us the set where the polynomials roots are
#These are the upper Cauchy's bound, the lower Cauchy's bound and the final set

function sidespan(enctor::AbstractVecor{<:Real})
    return 1+maximum(abs,enctor[begin:end-1])/enctor[end]
end
function midspan(enctor::AbstractVecor{<:Real})
    if enctor[1]==0
        return 0
    else
        return abs(enctor[1])/(abs(enctor[1])+maximum(abs,enctor[begin+1:end]))
    end
end
function rootspan(enctor::AbstractVector{<:Real})
    if enctor[1]==0
        println("The polynomials roots are in (", -sidespan(enctor), ",", sidespan(enctor), ")")
    else
        println("The polynomials roots are in (", -sidespan(enctor), ",", -midspan(enctor), ") U (", midspan(enctor), ",", sidespan(enctor), ")")
    end
end

#Now we can use use the previous functions to find usefull values that'll appear in Strum's Theorem
#We want to evaluate the signs of the elements of the Strum sequence at a,b, which are the sidespans
#and construct the sign vectors where the sign changes will be recorded

function strumstheorem(enctor::AbstractVector{<:Real})
    p0=Polynomials.Polynomial(enctor)
    p1=Polynomials.derivative(p0)
    signvector1=[sign(p0(-sidespan(enctor))), sign(p1(-sidespan(enctor)))]  #Here we record the signs at a 
    signvector2=[sign(p0(sidespan(enctor))), sign(p1(sidespan(enctor)))]    #Here we record the signs at b
    while length(p1.coeffs)>0
        p2=-rem(p0,p1)
        p0=p1
        p1=p2
        push!(signvector1, sign(p2(-sidespan(enctor))))
        push!(signvector2, sign(p2(sidespan(enctor))))
    end
    V1=0
    V2=0
    for i in 1:length(signvector1)-1
        if signvector1[i]!=signvector1[i+1]
            V1+=1
        end
        if signvector2[i]!=signvector2[i+1]
            V2+=1
        end
    end
    println("The number of the polynomal's roots is ", V1-V2)
end
function AddToJSON(enctor::AbstractVector{<:Real})
    ourdata=Dict("Polynomial"=>string(Polynomials.Polynomial(enctor)),
                 "Encoded vector"=>enctor, 
                 "Degree"=>length(enctor)-1, 
                 "Roots"=>roots(Polynomials.Polynomial(enctor))
                 )
    open("PolynomialDatabase.json", "w") do io
        JSON.print(io, ourdata, 4)                                           #CHANGE IF YOU ADD MORE INFORMATION TO THE DICTIONARY!!!!!!!!!!!!!!!!!!!
    end 
end
