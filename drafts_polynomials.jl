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
function parse_vec(s)
     if s== "Float64[]"
         return []
     end
     s=strip(s, ['[',']'])
     nums=split(s,",")
     return parse.(Float64,strip.(nums))
end 
  function descartes(coeffs)
      l=length(coeffs)
        pos_roots=0
       nonzero=[i for i in coeffs if i != 0]
        for i=1:l-1
            if nonzero[i]*nonzero[i+1]<0
                pos_roots=pos_roots+1
            end
        end
        return pos_roots
   end 