function f=plotBars(X)

fs=10; lw=2;
figure
colormap parula
hold on
bar(X','hist')

set(gca,'fontsize',fs)
xlabel('Intervention step (\tau)')
ylabel('X_i(\tau)')%,'rotation',0)
grid on
grid minor
box on