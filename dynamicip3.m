clear all
close all
clc

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end



[t,y]=ode15s(@ODE,[0,200],[(0.292241794328304+0.15)/2  0.018780544343014   0.011600697371167]);
hold on
plot(t,y,LineWidth=2)
legend("$c$ $(\mu M)$", '$n$ $(\%)$','$p$ $(\mu M)$')
xlabel('$t$')
set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,600,400])
set(gcf, 'Color', 'w')
hold off

function dydt=ODE(t,yy)
c=yy(1);
n=yy(2);
p=yy(3);

%parameters
Ve= 1;
Ke= 0.1;
Kflux= 4.89;
Kact= 0.2;
Hact= 2;
Hinh= 4;
HIP3= 4;
KIP3= 0.05;
Kinf= 52;
g= 0.5;
KPLC=0.2;
k3k=0.1;
K=0.3;
Vplc=1;
k5p=0.005;
alpha=0.15;
beta=4;
gamma=1/2;

%fluxes
Kinh= Kinf*(p^HIP3)/(p^HIP3+KIP3^HIP3);
PO1= ((beta*c-alpha)^Hact)/((beta*c-alpha)^Hact+Kact^Hact);
PO2= (Kinh^Hinh)/(Kinh^Hinh+(beta*c-alpha)^Hinh);
c1=(Ve*(beta*c-alpha)^2)/(Ke^2+(beta*c-alpha)^2);
eqp1 = (Vplc*(beta*c-alpha)^2)/(KPLC^2+(beta*c-alpha)^2);
eqp2 = (k3k*(beta*c-alpha)^2)/(K^2+(beta*c-alpha)^2);

dcdt=Kflux*n*PO1-c1;
dndt=g*(PO2-n);
dpdt=eqp1-(k5p+eqp2)*p;

dydt = gamma*[dcdt/beta;dndt;dpdt];
end