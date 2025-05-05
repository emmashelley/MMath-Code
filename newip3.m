clear all
close all
clc

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

[t,y]=ode15s(@ODE,[0,200],[(0.292241794328304+0.15)/2  0.9   0.011600697371167]);
hold on
plot(t,y,LineWidth=2)
legend("$c$ $(\mu M)$", '$n$ (Proportion)','$p$ $(\mu M)$')
xlabel('$t$ ($s$)')
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
Kinf= 2;
g= 0.5;
KPLC=0.2;
k3k=0.1;
K=0.4;
Vplc=1.3;
k5p=0.66;
alpha=0.15;
beta=4;
gamma=1/2;
K1=1;
m=4;
k=5;
g1=0.5;

%fluxes
Kinh= Kinf*(p^HIP3)/(p^HIP3+KIP3^HIP3);
PO1= ((beta*c-alpha)^Hact)/((beta*c-alpha)^Hact+Kact^Hact);
PO2= (Kinh^Hinh)/(Kinh^Hinh+(beta*c-alpha)^Hinh);
c1=(Ve*(beta*c-alpha)^2)/(Ke^2+(beta*c-alpha)^2);
eqp1 = (Vplc*(beta*c-alpha)^2)/(KPLC^2+(beta*c-alpha)^2);
eqp2 = (k3k*(beta*c-alpha)^2)/(K^2+(beta*c-alpha)^2);
PO3 = (K1^m)/(K1^m+p^m);

%ODEs
dcdt=Kflux*n*PO1*PO3-c1;
dndt=g*PO2-g1*n;
dpdt=eqp1-(k5p+eqp2)*p;

dydt = gamma*[dcdt/beta;dndt;dpdt];
end