function f=plotMultiOutb2(fb1,fb2,f1,f2,f3,fb1s,fb2s,f1s,f2s,f3s,tvec)

thresh=[12000,18000,24000];%[38172,53365,74966,82966];
hosp=[2000,12000,18000,24000,60000]+1000;%Points for text %change maxy accordingly
hospxloc=368;
%{
%Scenario A 6x1m:
deathsA=[.31123  1.8845  2.4248  3.4647  7.2388]*1e5;%Cumulative hospitalisations
gdpA=[6.6039 8.3682    8.5444    8.6105 8.8907]*1e5;
%Scenario B 6x1m:
deathsB=[1.1853  1.4516  2.1419  2.6059  7.0663]*1e5;%Cumulative hospitalisations
gdpB=[6.7032  7.1837    8.1432    8.1912   8.8907]*1e5;
%}
%
%Scenario A 3x2m:
deathsA=[.31123  1.8932    2.7590    3.5829  7.2388]*1e5;%Cumulative hospitalisations
gdpA=[6.6039 8.7331    8.7708    8.8021 8.8907]*1e5;
%Scenario B 3x2m:
deathsB=[1.1853  1.6620    2.2944    2.8687  7.0663]*1e5;%Cumulative hospitalisations
gdpB=[6.7032  8.3335    8.6337    8.7485   8.8907]*1e5;
%}
gdpAtext=cell(length(gdpA),1);
gdpBtext=gdpAtext;
for i=1:length(gdpAtext)
    gdpAtext{i}=strcat('£',num2str(round(gdpA(i)/1e3)),'bn');
    gdpBtext{i}=strcat('£',num2str(round(gdpB(i)/1e3)),'bn');
end

cmap=lines(7);
cmap2=.5*[1,1,1;1,1,1];%cmap+.5*(1-cmap);
figure
fs=10; lw=2;
maxY=max(fb2(:,3))+20000;
tld=find(fb1s(:,1)<tvec(3));
%
hold on
for i=1:length(thresh)
    plot([0,f1(end,1)],thresh(i)*[1,1],'-','linewidth',1,'color',0*[1,1,1])
end

for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
hh1=plot(fb1(:,1),fb1(:,3),'-','linewidth',lw,'color',cmap2(2,:));
hh2=plot(fb2(:,1),fb2(:,3),'--','linewidth',lw,'color',cmap2(2,:));
hh3=plot(f1(:,1),f1(:,3),'-','linewidth',lw,'color',cmap(1,:));%f1(:,1),f1(:,3),'-','linewidth',lw,'color',cmap(3,:));
hh4=plot(f2(:,1),f2(:,3),'linewidth',lw,'color',cmap(2,:));
hh5=plot(f3(:,1),f3(:,3),'linewidth',lw,'color',cmap(3,:));
plot(fb1s(tld,1),fb1s(tld,3),'k-','linewidth',lw);
%
lt=length(tvec);
points=tvec+10;
pointsy=.95*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',20);
text(tvec(2)+5,pointsy,'LD','fontsize',20);
for i=3:lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',20)
end

lh=length(hosp);
cmap3=[cmap(2,:);cmap(3:6,:);cmap(2,:)];
for i=1:lh
    text(hospxloc,hosp(i)+1000,gdpAtext{i})%cmap3(i,:));%230
end
xlabel('Time','FontSize',fs);
ylabel('Hospital occupancy','FontSize',fs);%yvar

xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})

set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])
%38172	53365   74966	82966
%18057   28115   38172	53365
legend([hh1,hh2,hh3,hh4,hh5],'LDA','FO','A(12,000)','A(18,000)','A(24,000)','location','west')
grid on
grid minor
box on
hold off
%%
figure
fs=10; lw=2;
maxY=max(fb2(:,3))+20000;
%
hold on
for i=1:length(thresh)
    plot([0,f1(end,1)],thresh(i)*[1,1],'-','linewidth',1,'color',0*[1,1,1])
end

for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
hh1=plot(fb1s(:,1),fb1s(:,3),'-','linewidth',lw,'color',cmap2(2,:));
hh2=plot(fb2s(:,1),fb2s(:,3),'--','linewidth',lw,'color',cmap2(2,:));
hh3=plot(f1s(:,1),f1s(:,3),'-','linewidth',lw,'color',cmap(1,:));
hh4=plot(f2s(:,1),f2s(:,3),'linewidth',lw,'color',cmap(2,:));
hh5=plot(f3s(:,1),f3s(:,3),'linewidth',lw,'color',cmap(3,:));
plot(fb1s(tld,1),fb1s(tld,3),'k-','linewidth',lw);
%
lt=length(tvec);
points=tvec+10;
pointsy=.95*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',20);
text(tvec(2)+5,pointsy,'LD','fontsize',20);
for i=3:lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',20)
end

lh=length(hosp);
cmap3=[cmap(2,:);cmap(3:6,:);cmap(2,:)];
for i=1:lh
    text(hospxloc,hosp(i)+1000,gdpBtext{i})
end
xlabel('Time','FontSize',fs);
ylabel('Hospital occupancy','FontSize',fs);%yvar

xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})

set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])
legend([hh1,hh2,hh3,hh4,hh5],'LDB','FO','B(12,000)','B(18,000)','B(24,000)','location','west')
grid on
grid minor
box on
hold off
%%
figure
gdpA=gdpA/1e3;
gdpB=gdpB/1e3;
hold on
plot([0,deathsA(end)-deathsA(1)],[0,gdpA(end)-gdpA(1)],':','linewidth',lw,'color',.5*[1,1,1])
h1=plot(deathsA(end)-deathsA(2:end),gdpA(end)-gdpA(2:end),'o-','color',.5*[0,1,0],'linewidth',lw','markersize',5,'markerfacecolor','w');%,.5*[0,1,0]);
plot(deathsA(end)-deathsA(1:2),gdpA(end)-gdpA(1:2),'o--','color',.5*[0,1,0],'linewidth',lw','markersize',5,'markerfacecolor','w');%,.5*[0,1,0]);
text(deathsA(end)-deathsA+5e3,gdpA(end)-gdpA-5,{'LDA','H_1','H_2','H_3','FO'})

plot([0,deathsB(end)-deathsB(1)],[0,gdpB(end)-gdpB(1)],':','linewidth',lw,'color',.5*[1,1,1])
h2=plot(deathsB(end)-deathsB(2:end),gdpB(end)-gdpB(2:end),'o-','color',.5*[1,0,0],'linewidth',lw','markersize',5,'markerfacecolor','w');%.5*[1,0,0]);
plot(deathsB(end)-deathsB(1:2),gdpB(end)-gdpB(1:2),'o--','linewidth',lw','color',.5*[1,0,0],'markersize',5,'markerfacecolor','w');%,.5*[1,0,0]);
%plot(deathsB(end)-deathsB(1),gdpB(end)-gdpB(1),'o','linewidth',lw,'markersize',15,'color',.5*[1,1,1])

text(deathsB(end)-deathsB(2:end)+5e3,gdpB(end)-gdpB(2:end)+5,{'H_1','H_2','H_3',''})
text(deathsB(end)-deathsB(1)+5e3,gdpB(end)-gdpB(1)+5,'LDB')
set(gca,'fontsize',fs)
xlabel('Hospitalisations averted');
ylabel('GDP loss (£bn)');
set(gca,'FontSize',fs);
legend([h1,h2],'Scenario A','Scenario B','location','NW')
axis tight%([0,deaths(end)-deaths(1),0,gdp(end)-gdp(1)])
grid on
grid minor
box on
hold off