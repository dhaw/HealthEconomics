function f=heCompareXoptimBars(x1,x2,x3,ddata,eecon)%(x1,x1s,x2,x2s,x3,x3s,x4,x4s)%x2 schools forced open
numInt=3;%Number of intervention periods
period=6/numInt;

plotTitle='Scenario A, 2-month periods';
NN=ddata.NNsector'/sum(ddata.NNsector);
legString={'Agriculture','Production','Construction','Distribution, transport, hotels and restaurants','Information and communication','Financial and insurance','Real estate','Professional and support activities','Government, health & education','Other services'};
%nvec=[0,cumsum(NN')];
x1s=max(1,ddata.xmin');%eecon.G\eecon.b/6);

%hmax=1;
H={'12,000','18,000','24,000'};

inds=repelem((1:10)',[3,23,1,9,4,3,1,9,4,6]');
inds=repmat(inds,1,3);

topn=5;
lx=length(x1)/numInt;
X=zeros(lx,3);
Y=zeros(topn,3);

x1=reshape(x1,lx,numInt);
%x1s=reshape(x1s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,1)=xdiffSum(:,2);
xdiffSum=sortrows(xdiffSum,2);
Y(:,1)=xdiffSum(1:topn,1);

x1=reshape(x2,lx,numInt);
%x1s=reshape(x2s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,2)=xdiffSum(:,2);
xdiffSum=sortrows(xdiffSum,2);
Y(:,2)=xdiffSum(1:topn,1);

x1=reshape(x3,lx,numInt);
%x1s=reshape(x3s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,3)=xdiffSum(:,2);
xdiffSum=sortrows(xdiffSum,2);
Y(:,3)=xdiffSum(1:topn,1);
%%
%GVA:
X=X.*eecon.obj*period;
%Employment:
%X=X.*repmat(NN,1,3);
%Y/Z not modified yet
%%
Z=Y;
Z(:,1)=X(Y(:,1),1);
Z(:,2)=X(Y(:,2),2);
Z(:,3)=X(Y(:,3),3);

%{
W=nan(size(X));
W(Y(:,1),1)=X(Y(:,1),1);
W(Y(:,2),2)=X(Y(:,2),2);
W(Y(:,3),3)=X(Y(:,3),3);
%}
f=Z;

%numInds=5;
%inds=xdiffSum(1:numInds,1);
%indSch=55;%find(xdiffSum(:,1)==55);


fs=10; lw=2; ms=8;
col1=.5*[1,1,1];
%cmap=.8*flipud(winter(numInds));
cmap=jet(10);
cmap=[0*[1,1,1];.4*[1,1,1];.8*[1,1,1];lines(7);];

maxy=500;%max(max(X));
miny=1000*floor(min(min(X))/1000);
figure
hold on
%{
plot([-1,6],[0,0],'k-','linewidth',1)
violinplot(X/6,H);
plot((1:4),Z/6,'kx','markersize',ms,'linewidth',lw,'markerfacecolor','w')
%}
%
%bar(X'/6,'facecolor',col1,'edgecolor',col1);%,'color','k');
%bar(W'/6);%,'color','k');
plot([1.5,1.5],[maxy,miny],'k--','linewidth',1.5)
plot([2.5,2.5],[maxy,miny],'k--','linewidth',1.5)
%plot([3.5,3.5],[maxy,miny],'k--','linewidth',1.5)
h=gobjects(1,10);
for i=1:10
    W=nan(size(X));
    W(inds==i)=X(inds==i);
    bar(W','facecolor',cmap(i,:),'edgecolor',cmap(i,:));
    h(i)=plot(-1,-1,'s','color',cmap(i,:),'markerfacecolor',cmap(i,:));
end
plot([.5,3.5],[0,0],'k-','linewidth',.1)
xticks(1:3)
xticklabels(H)
%}
%title(plotTitle)
xlabel('H_{max}')
ylabel('GVA loss (£m/6 months)')
set(gca,'fontsize',fs)
axis([.5,3.5,miny,maxy]) %500
legend(h,legString,'location','northeastoutside')

ax=get(gca);
ax.YAxis.Exponent=4;

grid on
grid minor
box on