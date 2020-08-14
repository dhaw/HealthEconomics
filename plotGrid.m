function f=plotGrid(x)
xvec=(1:8)';
xlab=num2str(xvec-2);
xvec=xvec-.5;
yvec=(1:10:63)';
ylab=num2str(yvec);
yvec=yvec-.5;

fs=10; lw=2;
figure
%colormap gray
imagesc(x)
xlabel('Period')
ylabel('Sector')
set(gca,'fontsize',fs,'xtick',xvec,'xticklabels',{'PRE','LD',xlab(3:end,:)},'ytick',yvec,'yticklabels',ylab)
axis([xvec(1),xvec(end)+1,.5,63.5])%yvec(1),yvec(end)+1])
caxis([0,1])
colorbar
grid on
box on
