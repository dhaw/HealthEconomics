function f=plotSensAnalRates(H,H2,f1)%(NNs64,ddatax,xoptim)
Hmax=38172;%38172	53365	74966	82966
fs=10; lw=2;
figure
hold on
plot([0,306],[Hmax,Hmax],'-','linewidth',lw,'color',.5*[1,1,1])
plot(1:306,f1(1:306,3),'k-','linewidth',lw)
%%
x=prctile(H,[5,50,95],2);
y1=x(:,1);
y2=x(:,2);
y3=x(:,3);
cmap=lines(7);
col1=cmap(2,:);
tvec=(1:306)';
tvec2=[tvec;flipud(tvec)];
inBetween=[y1;flipud(y3)];
maxy=max(y3);
fill(tvec2,inBetween,col1,'facealpha',.2);
plot(tvec,y1,'-','linewidth',1,'color',col1);
plot(tvec,y3,'-','linewidth',1,'color',col1);
h1=plot(tvec,y2,'-','linewidth',2,'color',col1);
%%
%
x=prctile(H2,[5,50,95],2);
y1=x(:,1);
y2=x(:,2);
y3=x(:,3);
cmap=lines(7);
col2=[0,0,0];%cmap(1,:);
tvec=(1:306)';
tvec2=[tvec;flipud(tvec)];
inBetween=[y1;flipud(y3)];
fill(tvec2,inBetween,col2,'facealpha',.2);
plot(tvec,y1,'--','linewidth',1,'color',col2);
plot(tvec,y3,'--','linewidth',1,'color',col2);
h2=plot(tvec,y2,'--','linewidth',2,'color',col2);
legend([h1,h2],'Full','Schools','location','NW')
%
%%
axis([0,306,0,maxy])
set(gca,'fontsize',fs)
xlabel('Time (days)')
ylabel('Hospital capacity')
grid on
grid minor
box on