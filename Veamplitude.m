clear all
close all
clc

VE=[0:0.01:1];
VPLC=[0:0.025:5];
n=length(VE);
m=length(VPLC);
amp=nan(n,m);
t_cutoff=200;

i=1;
tau=45;

for Ve=0:0.01:1
    j=1;
    for Vplc=0:0.025:5
        history = @(t) [(0.292241794328304+0.15)/2;  0.9;   0.011600697371167];
        sol=dde23(@(t,y,y_delayed)delay(t,y,y_delayed,Vplc,Ve),tau,history,[0 1000]);
        y1=sol.y(1,:);
        indices=sol.x>100;
        y_filtered=y1(indices);
        [y_max, idx]=max(y_filtered);
        amp(i,j)=y_max;
        j=j+1;
    end
    i=i+1;
end

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

hold on
x=linspace(0,1,size(amp,1));
y2=linspace(0,5,size(amp,2));
contourf(x,y2,amp','LineColor','none')
colormap('gray')
xlim([0 1])
ylim([0 5])
xlabel('$V_e$ ($\mu M/s$)')
ylabel('$V_{PLC}$ ($\mu M/s$)')
set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,500,375])
set(gcf, 'Color', 'w')
colorbar
ylabel(colorbar,'Maximum concentration of $c$ ($\mu M$)','Interpreter','latex')

load('vefreq.mat') 

osc=nan(n,m);

for a=1:n
    for b=1:m
        if frq(a,b)>0;
            osc(a,b)=1;
        else
            osc(a,b)=0;
        end
    end
end

highlight = double(osc);   
highlight(highlight == 0) = NaN;
hold on
contourf(x,y2,highlight',[1 1],'FaceColor','c','FaceAlpha',0.3,'LineColor','none')
hold on
contour(x,y2, highlight', [1 1], 'k', 'LineWidth', 2);
h = patch(NaN, NaN, 'c', 'FaceAlpha', 0.3);
legend(h, 'Oscillatory region');


function dydt=delay(t,y,y_delayed,Vplc,Ve)

c=y(1);
n=y(2);
p=y(3);

cdelay=y_delayed(1);
ndelay=y_delayed(2);
pdelay=y_delayed(3);

%parameters

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