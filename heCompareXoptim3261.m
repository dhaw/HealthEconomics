function f=heCompareXoptim3261(x1,x2,x3,s1,s2,s3,ddatax,eeconx)
obj=eeconx.obj;
x0=max(1,ddata.xmin');

x0=repmat(x0,1,6);
lx=length(x1)/6;
x1=(reshape(x1,lx,6)-1).*repmat(obj,1,6);
x1=sum(x1,1);
x2=(reshape(x2,lx,6)-1).*repmat(obj,1,6);
x2=sum(x2,1);
x3=(reshape(x3,lx,6)-1).*repmat(obj,1,6);
x3=sum(x3,1);
s1=(reshape(s1,lx,3)-1).*repmat(obj,1,3);
s1=reshape([s1;s1],lx,6);
s1=sum(s1,1);
s2=(reshape(s2,lx,3)-1).*repmat(obj,1,3);
s2=reshape([s2;s2],lx,6);
s2=sum(s2,1);
s3=(reshape(s3,lx,3)-1).*repmat(obj,1,3);
s3=reshape([s3;s3],lx,6);
s3=sum(s3,1);

fs=10; lw=2;
cmap=lines(7);%.8*flipud(winter(numInds));
figure
hold on
h=gobjects(6,1);%(numInds+1,1);
h(1)=plot(1:6,x1,'o-','linewidth',lw,'color',cmap(1,:),'markersize',5,'markerfacecolor','w');
h(2)=plot(1:6,s1,'o--','linewidth',lw,'color',cmap(1,:),'markersize',5,'markerfacecolor','w');
h(3)=plot(1:6,x2,'o-','linewidth',lw,'color',cmap(2,:),'markersize',5,'markerfacecolor','w');
h(4)=plot(1:6,s2,'o--','linewidth',lw,'color',cmap(2,:),'markersize',5,'markerfacecolor','w');
h(5)=plot(1:6,x3,'o-','linewidth',lw,'color',cmap(3,:),'markersize',5,'markerfacecolor','w');
h(6)=plot(1:6,s3,'o--','linewidth',lw,'color',cmap(3,:),'markersize',5,'markerfacecolor','w');
plot([1,6],[0,0],'k-','linewidth',1)
xlabel('Period')
ylabel('GVA loss (£m/6 months)')
set(gca,'fontsize',fs)
xticks(1:6)
miny=1000*floor(min(min([x1;s1]))/1000); maxy=500;
axis([1,6,miny,maxy])
lstr={'3x2m(12,000)','6x1m(12,000)','3x2m(18,000)','6x1m(18,000)','3x2m(24,000)','6x1m(24,000)'};
legend(h,lstr,'location','NE')%,'Hospitality')
grid on
grid minor
box on