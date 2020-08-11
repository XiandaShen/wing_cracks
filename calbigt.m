function[bigt]=calbigt(nc)
bigt=zeros(3,3,3,3);

for i=1:3
    for j=1:3
        for k=1:3
           for l=1:3
               djl=0;
               djk=0;
               dik=0;
               dil=0;
               
               if j==l
                   djl=1;
               end
               if j==k
                   djk=1;
               end
               if i==k
                   dik=1;
               end
               if i==l
                   dil=1;
               end
               
            bigt(i,j,k,l)=1/4*(nc(i)*nc(k)*djl+nc(i)*nc(l)*djk+nc(j)*nc(l)*dik+nc(j)*nc(k)*dil)-nc(i)*nc(j)*nc(k)*nc(l);
            
           end
        end
    end
end
end