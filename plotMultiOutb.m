function f=plotMultiOutb(fb1,fb2,f1,f2,f3,f4,tvec)
cmap=lines(7);
cmap2=cmap+.5*(1-cmap);

figure
fs=10; lw=2;
maxY=max(fb2(:,3))+20000;
%
hold on
thresh=[38172,53365,74966,82966];
for i=1:length(thresh)
    plot([0,f1(end,1)],thresh(i)*[1,1],'-','linewidth',lw','color',.5*[1,1,1])
end

for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
hh1=plot(fb1(:,1),fb1(:,3),'-','linewidth',lw,'color',cmap2(2,:));
hh2=plot(fb2(:,1),fb2(:,3),'--','linewidth',lw,'color',cmap2(2,:));
hh3=plot(f1(:,1),f1(:,3),'linewidth',lw,'color',cmap(3,:));
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
hosp=[1.3180e+04-5000,38172,53365,74966,82966,1.3074e+05-6000]+3000;%Points for text %change maxy accordingly
gdp={'753,590','788,730','834,060','863,200','870,420','1,035,047'};
lh=length(hosp);
cmap3=[cmap(2,:);cmap(3:6,:);cmap(2,:)];
for i=1:lh
    text(230,hosp(i),strcat('£',gdp(i),'m'),'color',cmap3(i,:));
end
xlabel('Time (days since 1st Jan)','FontSize',fs);
ylabel('Population','FontSize',fs);%yvar
set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])
%38172	53365   74966	82966
legend([hh1,hh2,hh3,hh4,hh5,hh6],'x=1','x=x_{min}','H_{max}=38172','H_{max}=53365','H_{max}=74966','H_{max}=82966','location','west')
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
grid on
grid minor
box on
hold off