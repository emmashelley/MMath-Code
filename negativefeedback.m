clear all
close all
clc

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

hold on
[t,y]=ode45(@ODE,[0, 200],[0.1, 0.7, 0.1]);
plot(t,y,'LineWidth',2)

legend('$c$ ($\mu M$)','$n$ (Proportion)','$p$ ($\mu M$)')
xlabel('$t$ ($s$)')
ylim([0 1.5])
set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,500,350])
set(gcf, 'Color', 'w')

function A=ODE(t,y)
c=y(1);
n=y(2);
p=y(3);

%parameters
K3K=0.4;
k3k=0.1;
k5p=0;
KPLC=0;
beta=0.185;
Ve=0.25;
Ke=0.1;
ctot=2;
v1=7.4;
k2=0.00148;
Kact=0.2;
Kinh=0.3;
KIP3=0.13;
tn=6.6;
VPLC=0.00045;

%fluxes
j1=c/(Kact+c);
j2=p/(KIP3+p);
j3=(ctot-c)/beta;
jchannel1=v1*(n^3)*(j1^3)*(j2^3)+k2;
jchannel=jchannel1*(j3-c);
j4=(c^2)/(Ke^2+c^2);
jpump=Ve*j4;

m1=(Kinh+c)/Kinh;

q1=(c^2)/(KPLC^2+c^2);
q2=(c^2)/(K3K^2+c^2);
vdeg1=VPLC*q1;
vdeg2=p*(k5p+k3k*q2);

%ODEs
dcdt=jchannel-jpump;
dndt=(1/tn)*(1-n*m1);
dpdt=vdeg1-vdeg2;

A=[dcdt; dndt; dpdt];
end