function f=plotGrid(x)
xvec=(1:6)';
xlab=num2str(xvec);
xvec=xvec-.5;
yvec=(1:10:63)';
ylab=num2str(yvec);
yvec=yvec-.5;

fs=10; lw=2;
figure
colormap gray
imagesc(x)
xlabel('Month')
ylabel('Sector')
set(gca,'fontsize',fs,'xtick',xvec,'xticklabels',xlab,'ytick',yvec,'yticklabels',ylab)
axis([xvec(1),xvec(end)+1,yvec(1),yvec(end)+1])
caxis([0,1])
grid on
box on
