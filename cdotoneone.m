function[r]=cdotoneone(a,b)
r=zeros(3,3);
for i=1:3
    for j=1:3
            r(i,j)=a(i)*b(j);
    end
end
end