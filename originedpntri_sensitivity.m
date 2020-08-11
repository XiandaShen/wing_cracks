mus=[0,0.2,0.4];
cstep=500;
cycle=1;
tstep=cycle*cstep*2; % time step
aepsilon=zeros(tstep,3);
lepsilon=zeros(tstep,3);
astress=zeros(tstep,3);
damage1=zeros(tstep,3);
damage2=zeros(tstep,3);
damage3=zeros(tstep,3);

damagem1=zeros(tstep,3);
damagem2=zeros(tstep,3);
damagem3=zeros(tstep,3);
damagew1=zeros(tstep,3);
damagew2=zeros(tstep,3);
damagew3=zeros(tstep,3);

for senn=1:3
% variable initiation
clearvars -except cstep cycle tstep senn mus aepsilon lepsilon astress damage1 damage2 damage3 damagem1 damagem2 damagem3 damagew1 damagew2 damagew3
Em=25000; %young modulus mpa
nu=0.32; % poisson ratio

Ec=18000; %young modulus mpa of crack


so=16/3*(1-nu^2)/Ec; % compliance of normal cracks
s1=32/3*(1-nu^2)/Ec/(2-nu)^2; % compliance of shear cracks


numc=37; % number of major cracks

Ko=80;  % constant for kic mode 1, mpa/l
sigmac=100; % in Kic
Ko2=344;  % constant for kic mode 2, mpa/l
sigmac2=319; % in Kic
co=4; % cohesion 
mu=mus(senn); % friction angel
cc1=400;
cc2=15;
gamap=0;

% cstep=zeros(cycle,1);
% cstep(1)=1000;
% cstep(2)=1050;
% cstep(3)=1100;

Nmc=zeros(numc,1); % normal coefficient
Tmc=zeros(numc,1); % shear coefficient
Bt=zeros(numc,1); % mainly use is friction reducion


am=zeros(numc,1); %length of major cracks
aw=zeros(numc,1); % length of wing cracks
ram=zeros(numc,tstep); % length of main cracks
raw=zeros(numc,tstep); % length of wing cracks

dam=zeros(numc,1); %length change of major cracks
daw=zeros(numc,1); % length change of wing cracks
Kicw=zeros(numc,1);


for i=1:numc
    am(i)=0.022;
    aw(i)=1e-10;
end

tvolume=1/100; % volume of rev/number of cracks


% dsigmam=[0 0 0; 0 0 0; 0 0 -0.004];
sigmamt=cell(cycle,1);
sigmamt{1}=[0 0 0; 0 0 0; 0 0 -32.269]; % global total stress
if senn==2
    sigmamt{1}=[0 0 0; 0 0 0; 0 0 -41.611]; % global total stress
end
if senn==3
        sigmamt{1}=[0 0 0; 0 0 0; 0 0 -51.151]; % global total stress
end
sigmamt{2}=[0 0 0; 0 0 0; 0 0 -36.5]; % global total stress
sigmamt{3}=[0 0 0; 0 0 0; 0 0 -39.28]; % global total stress
sigmamt{4}=[0 0 0; 0 0 0; 0 0 -39.62]; % global total stress
sigmamt{5}=[0 0 0; 0 0 0; 0 0 -39.79]; % global total stress
sigmamt{6}=[0 0 0; 0 0 0; 0 0 -39.928]; % global total stress
sigmamt{7}=[0 0 0; 0 0 0; 0 0 -40.0]; % global total stress
% sigmamt{1}=[0 0 0; 0 0 0; 0 0 -42.2]; % global total stress
%sigmamt{2}=[0 0 0; 0 0 0; 0 0 -42.01]; 
%sigmamt{3}=[0 0 0; 0 0 0; 0 0 -42.07]; 

sigmam=[-1 0 0; 0 -1 0; 0 0 -1]; % current global stress
sigmacn=zeros(numc,1);  % normal stress of cracks
sigmacs=zeros(numc,1);  % shear stress of cracks
shearnm=zeros(numc,1);  % final shear stress of cracks consider friction
sigmacml=zeros(numc,1);  % plane direction stress of cracks
sigmanw=zeros(numc,1); % normal stress of wing cracks
sigmacsvector=cell(numc,1);

epsilonm=zeros(3,3); % global strain initiation

tstrainedm=zeros(3,3); % global elastic damage strain of main cracks
tstrainedw=zeros(3,3); % global elastic damage strain of wing cracks

tstrained=cell(tstep,1); % global elastic damage strain
tstrainin=cell(tstep,1); % global irreversible strain
tstraine=cell(tstep,1); % global elastic strain
tstrain=cell(tstep,1); % global strain
tstress=cell(tstep,1); % global stress
for te=1:tstep
    tstrained(te)={zeros(3,3)};
    tstrainin(te)={zeros(3,3)};   
    tstraine(te)={zeros(3,3)};
    tstrain(te)={zeros(3,3)};
    tstress(te)={zeros(3,3)};   
end

% orientation of cracks
%nc=cell(numc,1);
nc=orientation37(numc);
mc=cell(numc,1);
mcmc=cell(numc,1);
ncnc=cell(numc,1);
for i=1:numc
    mc(i)={[0;0;0]};
end

% weight of orientation
wt=2*weight(numc);


rom=zeros(numc,1); % major cracks density
row=zeros(numc,1); % wing cracks density
drom=zeros(numc,1); % change of major cracks density
drow=zeros(numc,1); % change of wing cracks density

tdamage=cell(numc,1); % total damage 

for i=1:numc
    tdamage(i)={zeros(3,3)};
end

tdamagew=cell(numc,1); % total damage of wing crack

for i=1:numc
    tdamagew(i)={zeros(3,3)};
end

tdamagem=cell(numc,1); % total  damage of main crack

for i=1:numc
    tdamagem(i)={zeros(3,3)};
end




beta=zeros(numc,1); % normal displacement
gamma=cell(numc,1); % shear displacement

fdm=zeros(numc,1); % damage criterion for major cracks
fdm2=zeros(numc,1); % damage criterion for major cracks mode2
fdw=zeros(numc,1); % damage criterion for wing cracks

modec=zeros(numc,1);

tol=1e-12;  % newton criterion     
bign=cell(numc,1);
bigt=cell(numc,1);

rBt=zeros(tstep,i);
rshearnm=zeros(tstep,i);
rsigmacml=zeros(tstep,i);
rmodec=zeros(tstep,i);




%% loading step
for sn=1:tstep
    
%    determine the loading or unloading conditions
    
    if mod(floor((sn-1)/cstep),2)==1
        loads=-1;
    else
        loads=1;        
    end

    dsigmam=sigmamt{floor((sn-1)/cstep/2)+1}/cstep*loads;
    sigmam=sigmam+dsigmam;

% if (sn-1)<cstep(1)||((cstep(1)*2)<=(sn-1) && (sn-1)<(cstep(1)*2+cstep(2)))||((cstep(1)*2+cstep(2)*2)<=(sn-1) && (sn-1)<(cstep(1)*2+cstep(2)*2+cstep(3)))
% 
% sigmam=sigmam+dsigmam;
%     
% else
%   sigmam=sigmam-dsigmam;
% end      


% local stress
% normal stress on cracks
for i=1:numc
    sigmacn(i)=nc{i}'*sigmam*nc{i};
end

% shear stress on cracks
for i=1:numc
    sigmacsvector{i}=sigmam*nc{i}-nc{i}'*sigmam*nc{i}*nc{i};
    sigmacs(i)=norm(sigmacsvector{i});
    if sigmacs(i)==0
        sigmacs(i)=1e-20;
    end
end




for i=1:numc
    
% crack damage mode: 0 tensile, shear less than c, 1 compressive without
% shear strain no history, 2 tensile with shear strain, 3 compressive with
% shear strain, 4 compressive without shear strain with history 
    if sigmacn(i)>=0 
        if sigmacs(i)<co && modec(i)<2
        modec(i)=0;
        else
        modec(i)=2;
        end
    else
        if sigmacs(i)<co && modec(i)<2
            modec(i)=1;
        end
        if sigmacs(i)>=co && sigmacs(i)>=abs(mu*sigmacn(i))
            modec(i)=3;
        end
        if modec(i)>1 && sigmacs(i)>=abs(mu*sigmacn(i))
            modec(i)=3;
        end
        if sigmacs(i)>=co && sigmacs(i)<abs(mu*sigmacn(i))
            modec(i)=4;
        end        
        if sigmacs(i)<abs(mu*sigmacn(i)) && modec(i)>1
            modec(i)=4;
        end 
    end
    
% coefficient depends on mode
if modec(i)==0
    Nmc(i)=1;
    Tmc(i)=0;
    Bt(i)=1;
end
if modec(i)==1
    Nmc(i)=0;
    Tmc(i)=0;
    Bt(i)=1;
end    
if modec(i)==2
    Nmc(i)=1;
    Tmc(i)=1;
    Bt(i)=1;
end      
if modec(i)==3
    Nmc(i)=0;
    Tmc(i)=1;
    Bt(i)=(sigmacs(i)+mu*sigmacn(i))/sigmacs(i);
end  
if modec(i)==4
    Nmc(i)=0;
    Tmc(i)=1;
    Bt(i)=0;
end 


% Damage criterion
% major cracks
Kicm=am(i)^(3/2)/(1/Ko+am(i)/sigmac);
fdm(i)=sigmacn(i)*sqrt(pi*am(i))-Kicm;
oam=am(i);
if fdm(i)>0
am(i)=newam3(sigmacn(i),am(i),Ko,sigmac,tol);
end

% major cracks mode 2
Kicm2=am(i)^(3/2)/(1/Ko2+am(i)/sigmac2);
shearnm(i)=Bt(i)*sigmacs(i);
fdm2(i)=shearnm(i)*sqrt(pi*am(i))-Kicm2;
if fdm2(i)>0
am(i)=newam3(shearnm(i),am(i),Ko2,sigmac2,tol);
end



dam(i)=am(i)-oam;


% wing cracks

Kicw(i)=aw(i)^(3/2)/(1/Ko+aw(i)/sigmac);
mc(i)={sigmacsvector{i}/sigmacs(i)};
sigmacml(i)=mc{i}'*sigmam*mc{i};
shearnm(i)=Bt(i)*sigmacs(i);
fdw(i)=sqrt(pi)*(am(i)^(2)/aw(i)^(3/2)*shearnm(i)+sigmacml(i)*aw(i)^(1/2))-Kicw(i);
oaw=aw(i);


rshearnm(sn,i)=am(i)^(2)/aw(i)^(3/2)*shearnm(i);
rsigmacml(sn,i)=sigmacml(i)/shearnm(i);

sigmanw=0;
% sigmanw=zeros(3,3);    
% if sigmacn(i)<0  % wing crack does not develop when tension at main crack
    if fdw(i)>0.001
        aw(i)=newaw(sigmacml(i),shearnm(i),am(i),aw(i),Ko,sigmac,tol);
    end
    sigmanw=(am(i)/aw(i))^2*Bt(i)*sigmacs(i)+sigmacml(i);
% end
rBt(sn,i)=Bt(i);


rmodec(sn,i)=modec(i);

daw(i)=aw(i)-oaw;




% crack density
orom=rom(i);
orow=row(i);

rom(i)=am(i)^3/tvolume;
row(i)=aw(i)^3/tvolume;

drom(i)=rom(i)-orom;
drow(i)=row(i)-orow;

% normal displacement for major crack
beta(i)=rom(i)*so*Nmc(i)*sigmacn(i);

% shear displacement
gamma{i}=rom(i)*s1*Tmc(i)*Bt(i)*(sigmam*nc{i}-nc{i}'*sigmam*nc{i}*nc{i});

end





tstrainedm=zeros(3,3); % global elastic damage strain of main cracks
tstrainedw=zeros(3,3); % global elastic damage strain of wing cracks
dtstraininm=zeros(3,3); % global irreversible damage strain of main cracks
dtstraininw=zeros(3,3); % global irreversible damage strain of wing cracks
tddamage=zeros(3,3); % total damage increase at each step
tddamagem=zeros(3,3); % total main damage increase at each step
tddamagew=zeros(3,3); % total wing damage increase at each step

for i=1:numc

% to calculate main crack ed strain    
bign(i)={calbign(nc{i})};

bigt(i)={calbigt(nc{i})};

dsnmds=caldsnmds(nc{i}); % normal stress on main crack-----d sigma 2rd
dstmds=caldstmds(nc{i}); % shear stress on main crack-----d sigma 3rd
dnormstmds=1/sigmacs(i)*onedOT(sigmacsvector{i},dstmds); % norm of shear stress on main crack-----d sigma 2rd


    if mode(i)~=3
       dbtds=zeros(3,3);  % dBt----dsigma  2rd
    else
       dbtds=mu/sigmacs(i)*(dsnmds-1/(sigmacs(i))^2*sigmacn(i)*onedOT(sigmacsvector{i},dstmds));
    end

 tempstrainmc=rom(i)*wt(i)*(so*Nmc(i)*doubledotft(bign{i},sigmam)+s1*Tmc(i)*Bt(i)*doubledotft(bigt{i},sigmam)...
     +0.5*s1*dbtds*Tmc(i)*doubledottt(doubledottf(sigmam,bigt{i}),sigmam)); 
tstrainedm=tstrainedm+tempstrainmc;
 
% to calculate wing crack ed strain
dmds=1/sigmacs(i)*dstmds-1/sigmacs(i)^2*cdotonetwo(sigmacsvector{i},dnormstmds);   % 3rd
dmmds=cdotthreeone(dmds,mc{i})+cdotonethree(mc{i},dmds); % 4th
mtm=cdotoneone(mc{i},mc{i});    % 2rd
dsmlds=doubledotft(symidendityf,mtm)+doubledottf(sigmam,dmmds); % 2rd
dsnwds=(am(i)/aw(i))^2*(sigmacs(i)*dbtds+Bt(i)*dnormstmds)+dsmlds;  % 2rd
dsmmds=cdottwotwo(dsnwds,mtm)+sigmanw*dmmds;  % 4th
dssmmds=doubledotft(symidendityf,(sigmanw*mtm))+doubledottf(sigmam,dsmmds); %2rd

tempstrainwc=row(i)*wt(i)*so*(dssmmds-sigmanw*dsnwds);
if norm(tempstrainwc)>1
   tempstrainwc=zeros(3,3);
end
tstrainedw=tstrainedw+tempstrainwc;


% % to calculate main crack in strain 
% demcdrm=wt(i)*(so*Nmc(i)*doubledotft(bign{i},sigmam)+s1*Tmc(i)*Bt(i)*doubledotft(bigt{i},sigmam)...
%      +0.5*s1*dbtds*Tmc(i)*doubledottt(doubledottf(sigmam,bigt{i}),sigmam)); 
% 
% damdrm=tvolume/3*(rom(i))^(-2/3);
% dsnwdrm=2*am(i)/aw(i)^2*Bt(i)*sigmacs(i)*damdrm; % scalar
% dsnwdsdrm=2*am(i)/aw(i)^2*(sigmacs(i)*dbtds+Bt(i)*dnormstmds)*damdrm; % 2rd
% dsmmdrm=dsnwdrm*mtm; % 2rd
% dsmmdsdrm=cdottwotwo(dsnwdsdrm,mtm)+dsnwdrm*dmmds;  % 4th
% dssmmdsdrm=doubledotft(symidendityf,dsmmdrm)+doubledottf(sigmam,dsmmdsdrm); %2rd
% dewcdrm=wt(i)*row(i)*so*(dssmmdsdrm-dsnwdrm*dsnwds-sigmacs(i)*dsnwdsdrm); % 2rd
% 
% dedrm=demcdrm+dewcdrm;
% 
% dtstraininm=dtstraininm+cc1*dedrm*drom(i);
% 
% % to calculate wing crack in strain 
% dawdrw=tvolume/3*(row(i))^(-2/3);
% dsnwdrw=-2*am(i)^2/aw(i)^3*Bt(i)*sigmacs(i)*dawdrw;
% dsnwdsdrw=-2*am(i)^2/aw(i)^3*(sigmacs(i)*dbtds+Bt(i)*dnormstmds)*dawdrw; % 2rd
% dsmmdrw=dsnwdrw*mtm; % 2rd
% dsmmdsdrw=cdottwotwo(dsnwdsdrw,mtm)+dsnwdrw*dmmds;  % 4th
% dssmmdsdrw=doubledotft(symidendityf,dsmmdrw)+doubledottf(sigmam,dsmmdsdrw); %2rd
% dedrw=wt(i)*row(i)*so*(dssmmdsdrw-dsnwdrw*dsnwds-sigmacs(i)*dsnwdsdrw)...
%     +wt(i)*so*(dssmmds-sigmanw*dsnwds); % 2rd
% dtstraininw=dtstraininw+cc2*dedrw*drow(i);
% 
ddamage=drom(i)*cdotoneone(nc{i},nc{i})+drow(i)*cdotoneone(mc{i},mc{i});
ddamagem=drom(i)*cdotoneone(nc{i},nc{i});
ddamagew=drow(i)*cdotoneone(mc{i},mc{i});

ncnc(i)={cdotoneone(nc{i},nc{i})};
mcmc(i)={cdotoneone(mc{i},mc{i})};


tddamage=ddamage+tddamage;
tddamagew=ddamagew+tddamagew;
tddamagem=ddamagem+tddamagem;
end


% record a
raw(:,sn)=aw;
ram(:,sn)=am;

% record damage
if sn>1
tdamage(sn)={tdamage{sn-1}+tddamage};
else
tdamage(sn)={tddamage};   
end

if sn>1
tdamagem(sn)={tdamagem{sn-1}+tddamagem};
else
tdamagem(sn)={tddamagem};   
end

if sn>1
tdamagew(sn)={tdamagew{sn-1}+tddamagew};
else
tdamagew(sn)={tddamagew};   
end



% total strain

% calculate elastic damage strain
tstrained(sn)={tstrainedm+tstrainedw};


% calculate plastic strain
if sn>1
    [dtstrainin,ngamap]=plastic4n(trace(tdamage{sn}),trace(tddamage),tstrainin{sn-1},sigmam,dsigmam,gamap);
    tstrainin{sn}=tstrainin{sn-1}+dtstrainin;
else
    [dtstrainin,ngamap]=plastic4n(trace(tdamage{sn}),trace(tddamage),zeros(3,3),sigmam,dsigmam,gamap); 
    tstrainin{sn}=dtstrainin;
end
    
gamap=ngamap;
% tstrainin{sn}=zeros(3,3);

% calculate purely elastic strain
tstraine(sn)={(1+nu)/Em*sigmam-nu/Em*trace(sigmam)*eye(3)};

tstrain(sn)={tstrained{sn}+tstrainin{sn}+tstraine{sn}};

tstress(sn)={sigmam};







end




%% plot



for i=1:tstep
    aepsilon(i,senn)=tstrain{i}(3,3);
    lepsilon(i,senn)=tstrain{i}(1,1);    
    astress(i,senn)=tstress{i}(3,3)+1;    
    damage1(i,senn)=tdamage{i}(1,1);  
    damage2(i,senn)=tdamage{i}(2,2); 
    damage3(i,senn)=tdamage{i}(3,3); 
    damagem1(i,senn)=tdamagem{i}(1,1);  
    damagem2(i,senn)=tdamagem{i}(2,2); 
    damagem3(i,senn)=tdamagem{i}(3,3); 
    damagew1(i,senn)=tdamagew{i}(1,1);  
    damagew2(i,senn)=tdamagew{i}(2,2); 
    damagew3(i,senn)=tdamagew{i}(3,3); 
end
end

figure('Name','stress-strain','NumberTitle','off')
cruves{1}=plot(aepsilon,astress,'-k','LineWidth',1);
hold on
set(gcf, 'position', [200 200 400 300]);
set(gca,'FontSize',14)
xlabel('Axial strain','fontsize',14)
ylabel('Axial stress','fontsize',14)


load('/Users/caihejushi/Documents/academic/OneDrive - Georgia Institute of Technology/PHD/image_analysis/revcrack/expstressd.mat');
load('/Users/caihejushi/Documents/academic/OneDrive - Georgia Institute of Technology/PHD/image_analysis/revcrack/expstraina.mat');
load('/Users/caihejushi/Documents/academic/OneDrive - Georgia Institute of Technology/PHD/image_analysis/revcrack/expstrainl.mat');

lnodamage_stress=astress(:,2);
lnodamage_strain=lepsilon(:,2);
anodamage_strain=aepsilon(:,2);
    save('/Users/caihejushi/Documents/academic/OneDrive - Georgia Institute of Technology/PHD/image_analysis/revcrack/lnodamage_stress.mat','lnodamage_stress');    
    save('/Users/caihejushi/Documents/academic/OneDrive - Georgia Institute of Technology/PHD/image_analysis/revcrack/lnodamage_strain.mat','lnodamage_strain');  
    save('/Users/caihejushi/Documents/academic/OneDrive - Georgia Institute of Technology/PHD/image_analysis/revcrack/anodamage_strain.mat','anodamage_strain');  



figure('Name','stress-strain','NumberTitle','off')
cruves{1}=plot(-aepsilon(:,1),-astress(:,1),'-k','LineWidth',1);
hold on
cruves{2}=plot(-aepsilon(:,2),-astress(:,2),'-r','LineWidth',1);
hold on
cruves{3}=plot(-aepsilon(:,3),-astress(:,3),'-b','LineWidth',1);
hold on
cruves{4}=plot(-lepsilon(:,1),-astress(:,1),'--k','LineWidth',1);
hold on
cruves{5}=plot(-lepsilon(:,2),-astress(:,2),'--r','LineWidth',1);
hold on
cruves{6}=plot(-lepsilon(:,3),-astress(:,3),'--b','LineWidth',1);
hold on
axis([ -0.012,0.012,0, 60])
set(gcf, 'position', [200 200 700 350]);
set(gca,'FontSize',14)
xlabel('Lateral strain  <--->  Axial strain    ','fontsize',14)
ylabel('Differential stress (MPa)','fontsize',14)
legend([cruves{1},cruves{2},cruves{3},cruves{4},cruves{5},cruves{6}],'\mu=0 (Axial)',...
    '\mu=0.2 (Axial)','\mu=0.4 (Axial)','\mu=0 (Lateral)',...
    '\mu=0.2 (Lateral)','\mu=0.4 (Lateral)','Location','southwest')

figure('Name','stress-strain','NumberTitle','off')
cruves{1}=plot(-lepsilon(:,1),-astress(:,1),'-k','LineWidth',1);
hold on
cruves{2}=plot(-lepsilon(:,2),-astress(:,2),'-r','LineWidth',1);
hold on
cruves{3}=plot(-lepsilon(:,3),-astress(:,3),'-b','LineWidth',1);
hold on
%axis([-0.025,0,0, 50])
set(gcf, 'position', [200 200 400 300]);
set(gca,'FontSize',14)
% axis([-0.0035, 0,0, 50])

xlabel('Lateral strain','fontsize',14)
ylabel('Differential stress (MPa)','fontsize',14)
legend([cruves{1},cruves{2},cruves{3}],'\mu=0','\mu=0.2','\mu=0.4','Location','southeast')

% 
figure('Name','stress-damage3','NumberTitle','off')
cruves{1}=plot(-astress(:,1),damage3(:,1),'-k','LineWidth',1);
hold on
cruves{2}=plot(-astress(:,2),damage3(:,2),'-r','LineWidth',1);
hold on
cruves{3}=plot(-astress(:,3),damage3(:,3),'-b','LineWidth',1);
hold on
axis([0,60,0,0.8])
set(gcf, 'position', [200 200 400 300]);
set(gca,'FontSize',14)
xlabel('Differential stress (MPa)','fontsize',14)
ylabel('\Omega_3','fontsize',14)
legend([cruves{1},cruves{2},cruves{3}],'\mu=0','\mu=0.2','\mu=0.4','Location','southeast')


figure('Name','stress-damage3_m','NumberTitle','off')
cruves{1}=plot(-astress(:,1),damagem3(:,1),'-k','LineWidth',1);
hold on
cruves{2}=plot(-astress(:,2),damagem3(:,2),'-r','LineWidth',1);
hold on
cruves{3}=plot(-astress(:,3),damagem3(:,3),'-b','LineWidth',1);
hold on
set(gcf, 'position', [200 200 400 300]);
set(gca,'FontSize',14)
axis([0,60,0,0.5])
xlabel('Differential stress (MPa)','fontsize',14)
ylabel('\Omega_{m3}','fontsize',14)
legend([cruves{1},cruves{2},cruves{3}],'\mu=0','\mu=0.2','\mu=0.4','Location','southeast')

figure('Name','stress-damage3_m','NumberTitle','off')
cruves{1}=plot(-astress(:,1),damagem1(:,1),'-k','LineWidth',1);
hold on
cruves{2}=plot(-astress(:,2),damagem1(:,2),'-r','LineWidth',1);
hold on
cruves{3}=plot(-astress(:,3),damagem1(:,3),'-b','LineWidth',1);
hold on
set(gcf, 'position', [200 200 400 300]);
set(gca,'FontSize',14)
axis([0,60,0,0.5])
xlabel('Differential stress (MPa)','fontsize',14)
ylabel('\Omega_{m1}','fontsize',14)
legend([cruves{1},cruves{2},cruves{3}],'\mu=0','\mu=0.2','\mu=0.4','Location','southeast')



figure('Name','stress-damage3_w','NumberTitle','off')
cruves{1}=plot(-astress(:,1),damagew3(:,1),'-k','LineWidth',1);
hold on
cruves{2}=plot(-astress(:,2),damagew3(:,2),'-r','LineWidth',1);
hold on
cruves{3}=plot(-astress(:,3),damagew3(:,3),'-b','LineWidth',1);
hold on
set(gcf, 'position', [200 200 400 300]);
set(gca,'FontSize',14)
axis([0,60,0,0.5])
xlabel('Differential stress (MPa)','fontsize',14)
ylabel('\Omega_{w3}','fontsize',14)
legend([cruves{1},cruves{2},cruves{3}],'\mu=0','\mu=0.2','\mu=0.4','Location','southeast')
