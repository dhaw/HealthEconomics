function f=heCompareXoptimBarsDelta(x1,x2,x3,x4,ddata,eecon)%(x1,x1s,x2,x2s,x3,x3s,x4,x4s)%x2 schools forced open
numInt=3;%Number of intervention periods
period=6/numInt;

plotTitle='Scenario A, 2-month periods';
NN=ddata.NNsector'/sum(ddata.NNsector);
legString={'Agriculture','Production','Construction','Distribution, transport, hotels and restaurants','Information and communication','Financial and insurance','Real estate','Professional and support activities','Government, health & education','Other services'};
%nvec=[0,cumsum(NN')];
x1s=max(1,ddata.xmin');

%hmax=1;
delta={'0.61','0.62','0.63','0.64'};%,'0.7'};
ld=length(delta);

inds=repelem((1:10)',[3,23,1,9,4,3,1,9,4,6]');
inds=repmat(inds,1,ld);

topn=5;
lx=length(x1)/numInt;
X=zeros(lx,ld);
Y=zeros(topn,ld);

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

x1=reshape(x4,lx,numInt);
%x1s=reshape(x3s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,4)=xdiffSum(:,2);
xdiffSum=sortrows(xdiffSum,2);
Y(:,4)=xdiffSum(1:topn,1);
%{
x1=reshape(x5,lx,numInt);
%x1s=reshape(x3s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,5)=xdiffSum(:,2);
xdiffSum=sortrows(xdiffSum,2);
Y(:,5)=xdiffSum(1:topn,1);
%}
%%
%GVA:
X=X.*repmat(eecon.obj,1,ld)*period;
fs=10; lw=2; ms=8;
col1=.5*[1,1,1];
%cmap=.8*flipud(winter(numInds));
cmap=jet(10);
cmap=[0*[1,1,1];.4*[1,1,1];.8*[1,1,1];lines(7);];
maxy=500;%max(max(X));
miny=1000*floor(min(min(X))/1000);
figure
hold on
plot([1.5,1.5],[maxy,miny],'k--','linewidth',1.5)
plot([2.5,2.5],[maxy,miny],'k--','linewidth',1.5)
plot([3.5,3.5],[maxy,miny],'k--','linewidth',1.5)
%plot([4.5,4.5],[maxy,miny],'k--','linewidth',1.5)
h=gobjects(1,10);
for i=1:10
    W=nan(size(X));
    W(inds==i)=X(inds==i);
    bar(W','facecolor',cmap(i,:),'edgecolor',cmap(i,:));
    h(i)=plot(-1,-1,'s','color',cmap(i,:),'markerfacecolor',cmap(i,:));
end
plot([.5,5.5],[0,0],'k-','linewidth',1)
xticks(1:5)
xticklabels(delta)
%}
%title(plotTitle)
xlabel('NPI multiplier (\delta)')
ylabel('GVA loss (£m/6 months)')
set(gca,'fontsize',fs)
axis([.5,4.5,miny,maxy])
legend(h,legString,'location','northeastoutside')
grid on
grid minor
box on