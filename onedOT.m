function[r]=onedOT(a,b)
r=zeros(3,3);
for i=1:3
    for j=1:3
        
        for k=1:3
            r(i,j)=a(i)*b(i,j,k)+r(i,j);
        end
    end
end
end