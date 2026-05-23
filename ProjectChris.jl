using DataFrames, Random, LinearAlgebra
#Just a demo code, must be upgraded, works only for low degree polynomials

# First we set our parameters for solving our problem

size=101  #Degree of the polynomial +1
enctor=[]   #Encoded polynomial as a vector
maxcoef=2*rand()-1   #For Cauchy's bound, which is calculated later. We take advantage of Uni(0,1) for coeffs in [-1,1]
flip=0   #For the number of sign changes, calculatd for Descartes' rule of signs 
rc=0    #Root counter

#We now start constructing our vector that represents our polynomial. We'll use "coin" to determine which coefficients are 0
#The rest of the coefficients will be between -1 and 1

push!(enctor, maxcoef)   #Helps using a single loop in the next line
for i in 2:size-1
    coin=rand()
    if coin<0.5
       push!(enctor, 0)
    else push!(enctor, 2*rand()-1)
    end
    if abs(enctor[i])>=maxcoef
       global maxcoef=abs(enctor[i])
    end
    if enctor[i]*enctor[i-1]<0   #Sign changes counter
       global flip=flip+1
    end
end
For_highest_order_coefficient=rand()
while For_highest_order_coefficient==0.5
      global For_highest_order_coefficient=rand()
end
push!(enctor,2*For_highest_order_coefficient-1)   #Highest order coefficient
return enctor

#Next we have to determine a set where the roots of the polynomial are. We use Cauchy's bound for it with a small twist

sidespan=1+maxcoef/enctor[size]   #Upper Cauchy's bound
if enctor[1]==0
    midspan=0
else midspan=1/(1+max(maxcoef,enctor[size])/abs(enctor[1]))   #Lower Cauchy's bound   
end

#That way we have to partition the set (-sidespan, -midspan)U(midspan, sidespan), if 0 is not a root or 
#(-sidespan, sidespan) if 0 is a root. For accurate partition we'll use Mignotte's bound

formignotte=dot(enctor,enctor)   #"Polynomial's" norm_2 squared
mignotte=abs(enctor[size])*sqrt(3/(((size-1)^2)*(formignotte^(size-2))))   #Mignottes's bound. NOTE: THIS WAS CALCULATED USING COEFFS IN [-1,1]!!! CHANGE OTHERWIZE
loophelper=(sidespan-midspan)/mignotte   #this helps in the next loop
if enctor[1]==0    #That happens iff midspan=0. We don't need time consuming multiplications with 0s like in "else" case
    for i in 1:2*loophelper    
        b1=0   
        b2=0   
        for j in 1:size
            b1=b1+enctor[j]*(-sidespan+(i-1)*mignotte)^j   #f(x_1) for Bolzano's theorem
            b2=b2+enctor[j]*(-sidespan+i*mignotte)^j   #f(x_2) for Bolzano's theorem 
        end
        if b1*b2<0
            println("This polynomial has a root between ", -sidespan+(i-1)*mignotte, " and ", -sidespan+i*mignotte, " ")
            global rc=rc+1
        end
        if  rc==flip
            break
        end
    end
else for i in 1:2*loophelper    
        b1=0   #f(x_1) for Bolzano's theorem
        b2=0   #f(x_2) for Bolzano's theorem 
        for j in 1:size
            b1=b1+enctor[j]*(-sidespan+(i-1)*mignotte+floor(i/loophelper)*2*midspan)^j   #floor provides a "leap" if our root set is not connected in IR's usual topology
            b2=b2+enctor[j]*(-sidespan+i*mignotte+floor(i/loophelper)*2*midspan)^j   
        end
        if b1*b2<0
           println("This polynomial has a root between ", -sidespan+(i-1)*mignotte+floor(i/loophelper)*2*midspan, " and ", -sidespan+i*mignotte+floor(i/loophelper)*2*midspan, " ")
           global rc=rc+1
        end
        if  rc==flip
            break
        end
    end
end    
println("The polynomial has ", rc, " roots.")
