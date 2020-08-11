function[bign]=calbign(nc)
bign=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
           for l=1:3
            bign(i,j,k,l)=nc(i)*nc(j)*nc(k)*nc(l);
           end
        end
    end
end
end
