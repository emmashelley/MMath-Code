clear all
close all
clc

hold on

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

[t,y] = ode45(@politipos, [0,200],[1,0.2,0.5]);
plot(t,y(:,1),'-',LineWidth=2)
plot(t,y(:,2),'-',LineWidth=2)
plot(t,y(:,3),'-',LineWidth=2)
ylim([0 3.5])
legend('$c$ ($\mu M$)','$n$ (Proportion)','$p$ ($\mu M$)')
xlabel('$t$ ($s$)')
set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,500,350])
set(gcf, 'Color', 'w')
hold off

function A = politipos(t,y)
%parameters
Kact = 0.08;
KIP3 = 0.13;
v1 = 1.11;
Ve = 0.9;
Ke = 0.1;
KPLC = 0.2;
K3K = 0.4;
k5P = 0.66;
k3k = 0;
beta = 0.185;
tn = 12.5;
Kinh = 0.4;
k2 = 0.0203;
VPLC = 0.8;
ctot = 2;

%flux
Jchannel1 = y(1)/(Kact+y(1));
Jchannel2 = y(3)/(KIP3+y(3));
Jchannel3 = (ctot-y(1))/beta;
Jchannel = (v1*((y(2)*Jchannel1*Jchannel2)^3)+k2)*(Jchannel3-y(1));
Jpump = (Ve*y(1)^2)/(Ke^2+y(1)^2);
p1 = (y(1)^2)/(KPLC^2 + y(1)^2);
p2 = (y(1)^2)/(K3K^2 +y(1)^2);


%ODEs: c is y(1), n is y(2), p is y(3)
dcdt = Jchannel-Jpump;
dndt = (1/tn)*(1-y(2)*((Kinh+y(1)/Kinh)));
dpdt = VPLC*p1-(k5P+k3k*p2)*(y(3));

A = [dcdt; dndt; dpdt];
end
