function[dstmds]=caldstmds(nc)
dstmds=zeros(3,3,3);
for i=1:3
    for m=1:3
        for k=1:3
            dim=0;
            dik=0;
            if i==m
                dim=1;
            end
            if i==k
                dik=1;
            end
        dstmds(i,m,k)=0.5*(dim*nc(k)+dik*nc(m))-nc(i)*nc(m)*nc(k);
        end
    end
end
end