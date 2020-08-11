function[r]=cdotthreeone(a,b)
r=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
            r(i,j,k,l)=a(i,j,k)*b(l);
            end
        end
    end
end
end