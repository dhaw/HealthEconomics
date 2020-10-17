function f=heCompareXoptim(x1,x2,hmax,NN)%x2 schools forced open
%hmax=1;
titles={'H_{max}=18,057','H_{max}=28,115','H_{max}=38,172','H_{max}=53,365'};
lx=length(x1)/6;
x1=reshape(x1,lx,6);
x2=reshape(x2,lx,6);
xdiff=x2-x1;
xdiffSum=[(1:lx)',sum(xdiff,2)];

xdiffSum=sortrows(xdiffSum.*NN,2);
f=xdiffSum(1:20,:);

numInds=5;
inds=xdiffSum(1:numInds,1);
indSch=55;%find(xdiffSum(:,1)==55);

f=inds;

fs=10; lw=2;
cmap=.8*flipud(winter(numInds));
figure
hold on

xdiffx=xdiff;
xdiffx([inds;indSch],:)=[];
plot(xdiffx','-','linewidth',1,'color',.8*[1,1,1])
plot([0,10],[0,0],'k-','linewidth',1)

h=gobjects(numInds+1,1);
h(1)=plot(xdiff(indSch,:),'ko-','linewidth',2,'markersize',5,'markerfacecolor','w');%,'k');
for i=1:length(inds)
    h(i+1)=plot(xdiff(inds(i),:),'o-','linewidth',lw,'color',cmap(i,:),'markersize',5,'markerfacecolor','w');%,cmap(i,:));
end
%h(numInds+2)=plot(xdiff(56,:),'o--','linewidth',lw,'color',0*[1,1,1],'markersize',5,'markerfacecolor','w')%;,0*[1,1,1]);

xlabel('Month m')
ylabel('\Delta x_i(m)')
title(titles(hmax));
set(gca,'fontsize',fs)
xticks(1:6)
maxy=max(max(abs([xdiff(inds,:);xdiff(indSch,:)])));
axis([1,6,-.8,.5])%-maxy,maxy])
if hmax==1
    lstr={'Education','Sports/rec.','Arts & ent.','Air transp.','Paper','Furniture'};
elseif hmax==2
    lstr={'Education','Repair','Membership','Coke/petrol','Sports/rec.','Transport equip.'};
elseif hmax==3
    lstr={'Education','Travel agency','Sports/rec.','Textiles','Furniture','Basic metals'};
else
    lstr={'Education','Water coll.','Arts & ent.','Electric','Non-metallic','Sewerage'};
end
legend(h,lstr)%,'Hospitality')
grid on
grid minor
box on