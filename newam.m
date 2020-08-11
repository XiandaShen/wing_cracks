% use newton method
function[amn]=newam(sigmacn,am,Ko,sigmac,tol)
ttt=1;
while ttt==1
tfdm=sigmacn*sqrt(pi*am)-am^(3/2)/(1/Ko+am/sigmac);
dfdm=0.5*sigmacn*sqrt(pi)*am^(-0.5)-1.5*am^(0.5)/(1/Ko+am/sigmac)+am^(1.5)/sigmac/(1/Ko+am/sigmac)^2;

amn=am-tfdm/dfdm;

if norm(amn-am)<tol
    break
end
am=amn;


end
end