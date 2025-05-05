clear all
close all
clc


list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

[t,y] = ode45(@ODEs, [0,40], [4,1,20]);

plot(t,y(:,1),'-',LineWidth=2)
hold on
plot(t,y(:,2),'-',LineWidth=2)
plot(t,y(:,3),'-',LineWidth=2)
legend('$c$ ($\mu M$)','$n$ (Proportion)','$p$ ($\mu M$)')
xlabel('$t$ ($s$)')
set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,500,300])
set(gcf, 'Color', 'w')
ylim([0 35])
hold off


function A = ODEs(t,y)

%parameters
kflux = 4.8;
u0 = 0.567;
b = 0.111;
Ve = 20;
Ke = 0.06;
delta = 0;
a1 = 1;
a2 = 0.2;
Vp = 24;
tn = 2;
v4 = 6;
k4 = 1.1;
bosc = 0.08;
KIP3 = 4;
Kact = 0.7;
Kp = 0.4;
a = 0.97;
Kinh = 0.7;
ce = 14;
V1 = 0.889;

%fluxes: y(1) is c, y(2) is n, y(3) is p, 
Jchannel1 = (y(3)+u0*KIP3)/(KIP3+y(3));
Jchannel2 = (Kact*b+y(1))/(Kact+y(1));
Jchannel3 = ce-y(1);
Jchannel = kflux*Jchannel1*y(2)*Jchannel2*Jchannel3;

Jpump = (Ve*y(1))/(Ke+y(1));
Jin = a1+a2*v4;
Jpm = (Vp*y(1)^2)/(Kp^2+y(1)^2);

%ODEs
dcdt = Jchannel-Jpump+delta*(Jin-Jpm);
dpdt = v4*((y(1)+(1-a)*k4)/(y(1)+k4))-bosc*y(3);
dndt = (1/tn)*(((Kinh^2)/(Kinh^2+y(1)^2))-y(2));

A = [dcdt; dndt; dpdt];
end