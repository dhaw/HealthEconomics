function f=plotBars(X)
numPeriods=size(X,2);
fs=10; lw=2;
figure
colormap parula
hold on
bar(X','hist')
set(gca,'fontsize',fs)
xticks(1:numPeriods);
xticklabels({'PRE','LD','1','2','3','4','5','6'})
xlabel('Intervention step (\tau)')
ylabel('X_i(\tau)')%,'rotation',0)
grid on
grid minor
box on