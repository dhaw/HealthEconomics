function H=heSensAnalWfh(NNsx,ddatax,xoptim,tvec)
%
fact=(.8:.02:1);
lf=length(fact);
thresh=18000;
%
F=zeros(426,lf);
%G=F;
H=F;
parfor i=1:lf
    ddataxi=ddatax;
    ddataxi.wfhAv=fact(i)*ddataxi.wfhAv;
    [pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNsx,ddataxi);
    [f1,~]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,ddataxi);
    H(:,i)=f1(1:426,3);
end
cmap=flipud(jet(lf));
cmap2=.5*[1,1,1;1,1,1];%cmap+.5*(1-cmap);
figure
fs=10; lw=2;
maxY=35000;%max(max(H))+10000;
%tld=find(x<tvec(3));
%
hold on
plot([0,tvec(end)],thresh*[1,1],'-','linewidth',lw','color',.5*[1,1,1])
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for i=1:lf
    plot(1:tvec(end),H(:,i),'-','linewidth',lw,'color',cmap(i,:));
end
%plot(1:tvec(3),H(1:tvec(3),1),'k-','linewidth',lw);
%
lt=length(tvec);
points=tvec+10;
pointsy=.95*maxY;
txt={'1','2','3'};%,'4','5','6'};
text(10,pointsy,'PRE','fontsize',20);
text(tvec(2)+5,pointsy,'LD','fontsize',20);
for i=3:lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',20)
end
xlabel('Time','FontSize',fs);
ylabel('Hospital occupancy','FontSize',fs);%yvar

xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})

set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])
%legend([hh1,hh2,hh3,hh4,hh5],'LDA','U','A(12,000)','A(18,000)','A(24,000)','location','west')
grid on
grid minor
box on
hold off