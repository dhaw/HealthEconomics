function f=plotMultiOutb(fb1,fb2,f1,f2,f3,f4,fb1s,fb2s,f1s,f2s,f3s,f4s,tvec)
cmap=lines(7);
cmap2=cmap+.5*(1-cmap);
%%
figure
fs=10; lw=2;
maxY=max(fb2(:,3))+20000;
%
hold on
thresh=[18057,28115,38172,53365];%[38172,53365,74966,82966];
for i=1:length(thresh)
    plot([0,f1(end,1)],thresh(i)*[1,1],'-','linewidth',lw','color',.5*[1,1,1])
end

for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
hh1=plot(fb1(:,1),fb1(:,3),'-','linewidth',lw,'color',cmap2(2,:));
hh2=plot(fb2(:,1),fb2(:,3),'--','linewidth',lw,'color',cmap2(2,:));
hh3=plot(f1(:,1),f1(:,3),'-','linewidth',lw,'color',cmap(3,:));%f1(:,1),f1(:,3),'-','linewidth',lw,'color',cmap(3,:));
hh4=plot(f2(:,1),f2(:,3),'linewidth',lw,'color',cmap(4,:));
hh5=plot(f3(:,1),f3(:,3),'linewidth',lw,'color',cmap(5,:));
hh6=plot(f4(:,1),f4(:,3),'linewidth',lw,'color',cmap(6,:));
%
lt=length(tvec);
points=tvec+10;
pointsy=.95*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',20);
text(tvec(2)+5,pointsy,'LD','fontsize',20);
for i=3:8%lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',20)
end
hosp=[8000,18057,28115,38172,53365,110000]-4000;%Points for text %change maxy accordingly

%10 sector:
%gdp={'653,630','788,730','834,060','863,220','870,450','889,056'};
%gdp={'653,630','783,410','812,650','835,970','842,190','889,056'};%Schools opening

%64 sector
%Schools free:
gdp={'660,390',     '794,987',    '795,251',    '856,229',    '860,861',      '889,070'};

lh=length(hosp);
cmap3=[cmap(2,:);cmap(3:6,:);cmap(2,:)];
for i=1:lh
    text(246,hosp(i)+1000,strcat('£',gdp(i),'m'),'color','k')%cmap3(i,:));%230
end
xlabel('Time','FontSize',fs);
ylabel('Hospital occupancy','FontSize',fs);%yvar

xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})

set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])
%38172	53365   74966	82966
%18057   28115   38172	53365
legend([hh1,hh2,hh3,hh4,hh5,hh6],'x=x_{min}','x=1','H_{max}=18,057','H_{max}=28,115','H_{max}=38,172','H_{max}=53,365','location','west')
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
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
thresh=[18057,28115,38172,53365];%[38172,53365,74966,82966];
for i=1:length(thresh)
    plot([0,f1(end,1)],thresh(i)*[1,1],'-','linewidth',lw','color',.5*[1,1,1])
end

for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
hh1=plot(fb1s(:,1),fb1s(:,3),'-','linewidth',lw,'color',cmap2(2,:));
hh2=plot(fb2s(:,1),fb2s(:,3),'--','linewidth',lw,'color',cmap2(2,:));
hh3=plot(f1s(:,1),f1s(:,3),'-','linewidth',lw,'color',cmap(3,:));%f1(:,1),f1(:,3),'-','linewidth',lw,'color',cmap(3,:));
hh4=plot(f2s(:,1),f2s(:,3),'linewidth',lw,'color',cmap(4,:));
hh5=plot(f3s(:,1),f3s(:,3),'linewidth',lw,'color',cmap(5,:));
hh6=plot(f4s(:,1),f4s(:,3),'linewidth',lw,'color',cmap(6,:));
%
lt=length(tvec);
points=tvec+10;
pointsy=.95*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',20);
text(tvec(2)+5,pointsy,'LD','fontsize',20);
for i=3:8%lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',20)
end
hosp=[8000,18057,28115,38172,53365,110000]-4000;%Points for text %change maxy accordingly

%10 sector:
%gdp={'653,630','788,730','834,060','863,220','870,450','889,056'};
%gdp={'653,630','783,410','812,650','835,970','842,190','889,056'};%Schools opening

%64 sector
%Schools @80%+:
gdp={'670,320','712,065','730,533','760,733','819,030','889,070'};% //0356

lh=length(hosp);
cmap3=[cmap(2,:);cmap(3:6,:);cmap(2,:)];
for i=1:lh
    text(246,hosp(i)+1000,strcat('£',gdp(i),'m'),'color','k')%cmap3(i,:));%230
end
xlabel('Time','FontSize',fs);
ylabel('Hospital occupancy','FontSize',fs);%yvar

xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})

set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])
%38172	53365   74966	82966
%18057   28115   38172	53365
legend([hh1,hh2,hh3,hh4,hh5,hh6],'x=x_{min}','x=1','H_{max}=18,057','H_{max}=28,115','H_{max}=38,172','H_{max}=53,365','location','west')
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
grid on
grid minor
box on
hold off
%%
figure
%10 sector
%deaths=[.43489    .71695    .88321    1.1154    1.1848    1.5975]*1e5;
%gdp=[6.5363	7.8873    8.3406    8.6322    8.7045    8.89056]*1e5;

%64 sector
%Schools free:
%deaths=[1.0005  2.2178  2.7824  2.8911  3.6201  7.0663]*1e4;
deaths=[.12913  1.3742  1.9620  2.2543  3.0889  6.6682]*1e5;%Cumulative hospitalisations
gdp=[6.6039     7.9499    7.9525    8.5623    8.6086      8.8907]*1e5;
%Added point(s):
%deaths=[1.0005  2.2178 2.1507 2.7824 2.6711 2.8911  3.6201  7.0663]*1e4;
%gdp=[6.6039     7.9499  8.2422  7.9525  8.3347  8.5623    8.6086      8.8907]*1e5;
hold on
plot([0,deaths(end)-deaths(2)],[0,gdp(end)-gdp(2)],':','linewidth',lw,'color',.5*[1,1,1])
h1=plot(deaths(end)-deaths(2:end),gdp(end)-gdp(2:end),'o-','color',.5*[0,1,0],'linewidth',lw','markersize',5,'markerfacecolor','w');%,.5*[0,1,0]);
plot(deaths(end)-deaths(1:2),gdp(end)-gdp(1:2),'o--','color',.5*[0,1,0],'linewidth',lw','markersize',5,'markerfacecolor','w');%,.5*[0,1,0]);
%plot(deaths(end)-deaths(3),gdp(end)-gdp(3),'o','linewidth',lw,'markersize',15,'color',[1,0,0])
text(deaths(end)-deaths+1e3,gdp(end)-gdp,{'LD_A','H_1','H_2','H_3','H_4','Open'})
%text(deaths(end)-deaths(2:end)-3e3,gdp(end)-gdp(2:end)+5e3,{'H_1','H_1^+','H_2','H_3^-','H_3','H_4',''})

%Schools @80%+:
%deaths=[1.4768  2.1127  2.6089  3.1022  3.7323  7.0663]*1e4;
deaths=[.69034  1.3543  1.9048  2.4407  3.1716  7.0663]*1e5;%Cumulative hospitalisations
gdp=[6.7032     7.1207    7.3053    7.6073   8.1903      8.8907]*1e5;% // %6.7032
%hold on
plot([0,deaths(end)-deaths(2)],[0,gdp(end)-gdp(2)],':','linewidth',lw,'color',.5*[1,1,1])
h2=plot(deaths(end)-deaths(2:end),gdp(end)-gdp(2:end),'o-','color',.5*[1,0,0],'linewidth',lw','markersize',5,'markerfacecolor','w');%.5*[1,0,0]);
plot(deaths(end)-deaths(1:2),gdp(end)-gdp(1:2),'o--','linewidth',lw','color',.5*[1,0,0],'markersize',5,'markerfacecolor','w');%,.5*[1,0,0]);
plot(deaths(end)-deaths(1),gdp(end)-gdp(1),'o','linewidth',lw,'markersize',15,'color',.5*[1,1,1])
%plot(deaths(end)-deaths(2),gdp(end)-gdp(2),'x','linewidth',lw,'markersize',15,'color',[1,0,0])
%plot(deaths(end)-deaths(3),gdp(end)-gdp(3),'o','linewidth',lw,'markersize',15,'color',[1,0,0])

text(deaths(end)-deaths(2:end)-3e3,gdp(end)-gdp(2:end)+5e3,{'H_1','H_2','H_3','H_4',''})

text(deaths(end)-deaths(1)-4.5e3,gdp(end)-gdp(1),'LD_B')
set(gca,'fontsize',fs)
xlabel('Hospitalisations averted');
ylabel('GDP loss (£m)');
set(gca,'FontSize',fs);
legend([h1,h2],'Scenario A','Scenario B','location','NW')
axis tight%([0,deaths(end)-deaths(1),0,gdp(end)-gdp(1)])
grid on
grid minor
box on
hold off