clear all
close all
clc

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

tau=45;
history = @(t) [(0.092241794328304+0.15)/2;  0.9;   0.011600697371167];
sol=dde23(@delay,tau,history,[0 1000]);
hold on
plot(sol.x, sol.y,LineWidth=2)
legend("$c$ $(\mu M)$", '$n$ (Proportion)','$p$ $(\mu M)$','NumColumns',2)

set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,500,350])
set(gcf, 'Color', 'w')
xlabel('$t$ ($s$)')
ylim([0 4])
hold off

function dydt=delay(t,y,y_delayed)

c=y(1);
n=y(2);
p=y(3);

pdelay=y_delayed(3);

%parameters
Ve= 0.4;
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
K=0.3;
Vplc=1.7;
k5p=0.66;
alpha=0.15;
beta=4;
gamma=1/4;
K1=1;
m=4;
k=5;
g1=0.5;

%fluxes
Kinh= Kinf*(pdelay^HIP3)/(pdelay^HIP3+KIP3^HIP3);
PO1= ((beta*c-alpha)^Hact)/((beta*c-alpha)^Hact+Kact^Hact);
PO2= (Kinh^Hinh)/(Kinh^Hinh+(beta*c-alpha)^Hinh);
c1=(Ve*(beta*c-alpha)^2)/(Ke^2+(beta*c-alpha)^2);
eqp1 = (Vplc*(beta*c-alpha)^2)/(KPLC^2+(beta*c-alpha)^2);
eqp2 = (k3k*(beta*c-alpha)^2)/(K^2+(beta*c-alpha)^2);
PO3 = (K1^m)/(K1^m+pdelay^m);

%DDEs
dcdt=Kflux*n*PO1*PO3-c1;
dndt=g*PO2-g1*n;
dpdt=eqp1-(k5p+eqp2)*p;

dydt = gamma*[dcdt/beta;dndt;dpdt];
end