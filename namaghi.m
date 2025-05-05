clear all
close all
clc

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

[t,y]=ode45(@ODE,[0,100],[0.4,0.5]);
hold on
plot(t,y,LineWidth=2)
legend('$c$ ($\mu M$)','$n$ ($\%$)')
xlabel('$t$')
ylim([0 2])
set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,500,350])
set(gcf, 'Color', 'w')
hold off

function dydt=ODE(t,yy)

c=yy(1);
n=yy(2);

%parameters
Ve= 1;
Ke= 0.1;
Kflux= 4.89;
Kact= 0.2;
Hact= 2;
Hinh= 4;
HIP3= 4;
KIP3= 0.05;
p= 0.0085;
Kinf= 52;
g= 0.5;

%fluxes
Kinh= Kinf*(p^HIP3)/(p^HIP3+KIP3^HIP3);
PO1= ((10*c)^Hact)/((10*c)^Hact+Kact^Hact);
PO2= (Kinh^Hinh)/(Kinh^Hinh+(10*c)^Hinh);
c1=(Ve*(10*c)^2)/(Ke^2+(10*c)^2);

%ODEs
dcdt=Kflux*n*PO1-c1;
dndt=g*(PO2-n);

dydt = [dcdt/10;dndt];
end