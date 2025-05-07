clear all
close all
clc

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end


K=1;
n=4;
x=linspace(0,5);

y=(K^n)./(K^n+x.^n);

plot(x,y,'LineWidth',2,'DisplayName','$P_{O3}=\frac{K^4}{K^4+p^4}$')
hold on
x_point=1;
y_point=1/2;
scatter(x_point, y_point, 100, 'r', 'filled','DisplayName','$K=$ Half inactivation constant');
set(gca,'FontSize',24)
set(gcf, 'Position', [50,50,600,400])
set(gcf, 'Color', 'w')
xlabel('$p$ ($\mu M$)' )
ylabel('$P_{O3}$')
legend

