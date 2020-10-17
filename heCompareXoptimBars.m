function f=heCompareXoptimBars(x1,x2,x3,x4,ddata)%(x1,x1s,x2,x2s,x3,x3s,x4,x4s)%x2 schools forced open
NN=ddata.NNsector'/sum(ddata.NNsector);
legString={'Agriculture','Production','Construction','Distribution, transport, hotels and restaurants','Information and communication','Financial and insurance','Real estate','Professional and support activities','Government, health & education','Other services'};
%nvec=[0,cumsum(NN')];
x1s=1;
x2s=1;
x3s=1;
x4s=1;
%hmax=1;
H={'18,057','28,115','38,172','53,365'};

inds=repelem((1:10)',[3,23,1,9,4,3,1,9,4,6]');
inds=repmat(inds,1,4);

topn=5;
lx=length(x1)/6;
X=zeros(lx,4);
Y=zeros(topn,4);

x1=reshape(x1,lx,6);
%x1s=reshape(x1s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,1)=xdiffSum(:,2).*NN;%****
xdiffSum=sortrows(xdiffSum,2);
Y(:,1)=xdiffSum(1:topn,1);

x1=reshape(x2,lx,6);
%x1s=reshape(x2s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,2)=xdiffSum(:,2).*NN;%****
xdiffSum=sortrows(xdiffSum,2);
Y(:,2)=xdiffSum(1:topn,1);

x1=reshape(x3,lx,6);
%x1s=reshape(x3s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,3)=xdiffSum(:,2).*NN;%****
xdiffSum=sortrows(xdiffSum,2);
Y(:,3)=xdiffSum(1:topn,1);

x1=reshape(x4,lx,6);
%x1s=reshape(x4s,lx,6);
xdiff=x1-x1s;
xdiffSum=[(1:lx)',sum(xdiff,2)];
X(:,4)=xdiffSum(:,2).*NN;%****
xdiffSum=sortrows(xdiffSum,2);
Y(:,4)=xdiffSum(1:topn,1);

Z=Y;
Z(:,1)=X(Y(:,1),1);
Z(:,2)=X(Y(:,2),2);
Z(:,3)=X(Y(:,3),3);
Z(:,4)=X(Y(:,4),4);
%{
W=nan(size(X));
W(Y(:,1),1)=X(Y(:,1),1);
W(Y(:,2),2)=X(Y(:,2),2);
W(Y(:,3),3)=X(Y(:,3),3);
W(Y(:,4),4)=X(Y(:,4),4);
%}
f=Z;

%numInds=5;
%inds=xdiffSum(1:numInds,1);
%indSch=55;%find(xdiffSum(:,1)==55);


fs=10; lw=2; ms=8;
col1=.5*[1,1,1];
%cmap=.8*flipud(winter(numInds));
cmap=jet(10);
maxy=max(max(X))/6;
miny=min(min(X))/6;
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
plot([3.5,3.5],[maxy,miny],'k--','linewidth',1.5)
h=gobjects(1,10);
for i=1:10
    W=nan(size(X));
    W(inds==i)=X(inds==i);
    bar(W'/6,'facecolor',cmap(i,:),'edgecolor',cmap(i,:));
    h(i)=plot(-1,-1,'s','color',cmap(i,:),'markerfacecolor',cmap(i,:));
end
%plot([-1,6],[0,0],'k-','linewidth',1)
xticks(1:4)
xticklabels(H)
%}
xlabel('H_{max}')
ylabel('Closure (employment)')
set(gca,'fontsize',fs)
axis([.5,4.5,miny,maxy])
legend(h,legString,'location','northeastoutside')
grid on
grid minor
box on