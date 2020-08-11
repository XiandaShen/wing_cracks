clear;
c=4;
q=0:1:40;
alpha=20:2:100;
[qq,aa] = meshgrid(q,alpha);
p=-qq.^2./aa+c;

figure('Name','yield surface','NumberTitle','off')
surf(p,aa,qq)
set(gca,'FontSize',16)
xlabel('p (MPa)','fontsize',16)
ylabel('\alpha_p (MPa)','fontsize',16)
zlabel('q (MPa)','fontsize',16)