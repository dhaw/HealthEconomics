function f=heCompareXoptim(x1,hmax,ddatax,eeconx)%x2 schools forced open
NN=ddatax.NNsector;
obj=eeconx.obj;
monthPeriod=2;
numInt=6/monthPeriod;%3;
miny=-6000; maxy=100;
%hmax=1;
titles={'H_{max}=12,000','H_{max}=18,000','H_{max}=24,000'};
lx=length(x1)/numInt;
xdiff=(reshape(x1,lx,numInt)-1).*repmat(obj,1,numInt).*monthPeriod;
%
%Ranked:
xdiffSum=[(1:lx)',sum(xdiff,2)];
xdiffSum=sortrows(xdiffSum,2);
%f=xdiffSum(1:20,:);
numInds=6;
inds=xdiffSum(1:numInds,1);

inds(inds==55)=[];

f=inds;%Top 5
indSch=55;%find(xdiffSum(:,1)==55);
%
%Override with specific sectors:
inds=[27,30,31,32,33,36,58,59]';
numInds=length(inds);
%}
fs=10; lw=2;
cmap=lines(numInds);%.8*flipud(winter(numInds));
figure
hold on

xdiffx=xdiff;
xdiffx([inds;indSch],:)=[];
plot(xdiffx','-','linewidth',1,'color',.8*[1,1,1])
plot([0,numInt+1],[0,0],'k-','linewidth',1)

h=gobjects(5,1);%(numInds+1,1);
h(1)=plot(xdiff(indSch,:),'ko--','linewidth',2,'markersize',5,'markerfacecolor','w');%,'k');
h(2)=plot(sum(xdiff(inds(3:5),:),1),'o-','linewidth',lw,'color',cmap(1,:),'markersize',5,'markerfacecolor','w');
inds(3:5)=[];
for i=1:length(inds)%-3
    h(i+2)=plot(xdiff(inds(i),:),'o-','linewidth',lw,'color',cmap(i+1,:),'markersize',5,'markerfacecolor','w');%,cmap(i,:));
end
%h(numInds+2)=plot(xdiff(56,:),'o--','linewidth',lw,'color',0*[1,1,1],'markersize',5,'markerfacecolor','w')%;,0*[1,1,1]);

xlabel('Period m')
ylabel('GVA)')
%title(titles(hmax));
set(gca,'fontsize',fs)
xticks(1:numInt)
%maxy=max(max(abs([xdiff(inds,:);xdiff(indSch,:)])));
miny=-8000; maxy=100;
axis([1,numInt,miny,maxy])
lstr={'Education','Construction','Retail','Transport','Accom.','Arts & ent.s','Sports & rec.'};
%{
if hmax==1
    lstr={'Education','Sports/rec.','Arts & ent.','Air transp.','Paper','Furniture'};
elseif hmax==2
    lstr={'Education','Repair','Membership','Coke/petrol','Sports/rec.','Transport equip.'};
elseif hmax==3
    lstr={'Education','Travel agency','Sports/rec.','Textiles','Furniture','Basic metals'};
else
    lstr={'Education','Water coll.','Arts & ent.','Electric','Non-metallic','Sewerage'};
end
%}
legend(h,lstr,'location','SW')%,'Hospitality')
grid on
grid minor
box on