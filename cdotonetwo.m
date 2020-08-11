function[r]=cdotonetwo(a,b)
r=zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            r(i,j,k)=a(i)*b(j,k);
        end
    end
end
end
