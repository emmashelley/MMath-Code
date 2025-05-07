clear all
close all
clc

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

VPLC=[5:-0.5:0];
KFLUX=[2:2:10];
n=length(VPLC);
m=length(KFLUX);
tau=45;

t = tiledlayout(n, m);
t.TileSpacing = 'compact';
t.Padding = 'compact';
row=0;

for Vplc = VPLC
    row=row+1;
    col=0;
    for Kflux = KFLUX
        col=col+1;
        history = @(t) [(0.292241794328304+0.15)/2; 0.9; 0.011600697371167];
        sol = dde23(@(t,y,y_delayed) delay(t,y,y_delayed,Vplc,Kflux), tau, history, [0 1000]);
        nexttile
        plot(sol.x, sol.y,'LineWidth',1)
        xlim([0 1000])
        ylim([0 7])
        if col==1
            ylabel(sprintf('$V_{PLC} = %.1f$',Vplc),"FontSize",12,"FontWeight","normal")
        end
        if row==n
            xlabel(sprintf('$K_{flux} = %.1f$',Kflux),"FontSize",12,"FontWeight","bold")
        end
        if Vplc==5 && Kflux==2
            legend({'$c$','$n$','$p$'})
        end
    end
end


set(gcf, 'Position', [100, 50, 800, 1500])
set(gcf, 'Color', 'w')
export_fig kfluxsimulations.png -r600



function dydt=delay(t,y,y_delayed,Vplc,Kflux)

c=y(1);
n=y(2);
p=y(3);

cdelay=y_delayed(1);
ndelay=y_delayed(2);
pdelay=y_delayed(3);

%parameters
Ve= 0.4;
Ke= 0.1;
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


dcdt=Kflux*n*PO1*PO3-c1;
dndt=g*PO2-g1*n;
dpdt=eqp1-(k5p+eqp2)*p;

dydt = gamma*[dcdt/beta;dndt;dpdt];
end