function[dsnmds]=caldsnmds(nc)
dsnmds=zeros(3,3);
for i=1:3
    for j=1:3
        dsnmds(i,j)=nc(i)*nc(j);
    end
end
end
