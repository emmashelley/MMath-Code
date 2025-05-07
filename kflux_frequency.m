clear all
close all
clc


KFLUX=[0:0.1:10];
VPLC=[0:0.025:5];
n=length(KFLUX);
m=length(VPLC);
frq=nan(n,m);
t_cutoff=200;

i=1;

tau=45;
for Kflux=0:0.1:10
    j=1;
    for Vplc=0:0.025:5
        history = @(t) [(0.292241794328304+0.15)/2;  0.9;   0.011600697371167];
        sol=dde23(@(t,y,y_delayed)delay(t,y,y_delayed,Vplc,Kflux),tau,history,[0 2000]);
        t_start=500;
        indices=sol.x>t_start;
        signal=sol.y(1,indices)-mean(sol.y(1,indices));
        N=length(sol.x(indices));
        dt=mean(diff(sol.x(indices)));
        Fs=1/dt;
        f=(-N/2:N/2-1)*(Fs/N);
        power=abs(fftshift(fft(signal))).^2;
        [~, idx_max] = max(power);
        hz_freq = abs(f(idx_max));
        timescale = 1/hz_freq;
        frq(i,j)=hz_freq;
        j=j+1;
    end
    i=i+1;
end

frq(frq<0.002)=NaN;
frq=1./frq;
frq(frq<100)=NaN;
frq=1./frq;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

x=linspace(0,10,size(frq,1));
y2=linspace(0,5,size(frq,2));
contourf(x,y2,frq','LineColor','none')
xlim([0 10])
ylim([0 5])
xlabel('$K_{flux}$')
ylabel('$V_{PLC}$ ($\mu M/s$)')
set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,500,375])
set(gcf, 'Color', 'w')
colorbar
ylabel(colorbar,'Frequency of $c$ ($H_z$)','Interpreter','latex')
export_fig freq_kflux_final2.png -r600

save('kfluxfreq.mat','frq')

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