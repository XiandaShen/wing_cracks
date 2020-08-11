function[dep,gamapn]=plastic4n(d,dd,ep,sigma,dsigma,gamap)
%%  parameters
% input
% dd=0;
% ep=0;
% sigma=0;
% plastic
b=0;
B=0.00307;
C=4;
gx=0.5;
apo=20.25;
apm=490;
eta=-0.79;
dep=zeros(3,3);


% B=0.01;
% C=4;
% gx=0.5;
% apo=20.25;
% apm=490;
% eta=-0.5;
% dep=zeros(3,3);
%%
p=1/3*trace(sigma);
S=sigma-p*eye(3);
J2=1/2*doubledottt(S,S);
J3=det(S);
sintheta=-3*sqrt(3)/2*J3/(J2^(3/2));
gamapn=gamap;
dgamap=0;
q=sqrt(3*J2);
alphap=(1-gx*d)*(apo+(apm-apo)*gamap/(gamap+B));
h=1-b*sintheta;
f=q^2*h^2+alphap*(p-C);

%
if f>0
    alphap=-q^2/(p-C)
    gamapn=B/((apm-apo)/(alphap/(1-gx*d)-apo)-1);
    dgamap=gamapn-gamap;
    
    dqdj2=sqrt(3)/2*J2^(-0.5);
    dj2dsigma=S;
    dqdsigma=dqdj2*dj2dsigma;
    dpdsigma=1/3*eye(3);
    dgdsigma=dqdsigma-eta*(1-gx*d)*dpdsigma;
    dgamadlamda=sqrt(2/3*(doubledottt(dgdsigma,dgdsigma)-1/3*(trace(dgdsigma))^2));
    
    dlamda=dgamap/dgamadlamda;
    dep=dlamda*dgdsigma;
end
% if dgamap<0
%     dgamap=0;
% end

%f=q^2*h^2+alphap*(p-C)
