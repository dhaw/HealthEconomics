function f=plotMultiOut(f1,f2,tvec)
dodiff=1;
yvar='Inc./hosp.';%'Susceptibles'; 'Hospitalisations';
solvetype=2;
tend=tvec(end);%720;%For plot only
cmap=lines(7);
cmap2=cmap+.5*(1-cmap);
col2=.8*[1,1,1];

t1=f1(:,1); inc1=f1(:,2); h1=f1(:,3);%heSimCovid19 output
t2=f2(:,1); inc2=f2(:,2); h2=f2(:,3);
if dodiff==1
    inc1=-diff(inc1,1);%f1=Sout
    tdiff=diff(t1,1);
    inc1=inc1./tdiff;%repmat - older version?
    t1(1)=[];
    h1(1)=[];
    inc2=-diff(inc2,1);%f1=Sout
    tdiff=diff(t2,1);
    inc2=inc2./tdiff;%repmat - older version?
    t2(1)=[];
    h2(1)=[];
end
figure
fs=10; lw=2;
maxY=max([inc1;inc2;h1;h2]);
%
%Unlogged plots:
h=zeros(1,4);
hold on
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
hh3=plot(t2,inc2,'linewidth',lw,'color',cmap2(1,:));
hh4=plot(t2,h2,'--','linewidth',lw,'color',cmap2(2,:));
hh1=plot(t1,inc1,'linewidth',lw,'color',cmap(1,:));
hh2=plot(t1,h1,'--','linewidth',lw,'color',cmap(2,:));
%}
%{
%Logged plots:
hold on
semilogy(tout,Yall);
%}
lt=length(tvec);
points=tvec+10;
pointsy=.93*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',20);
text(tvec(2)+5,pointsy,'LD','fontsize',20);
for i=3:8%lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',20)
end

xlabel('Time','FontSize',fs);
ylabel('Population','FontSize',fs);%yvar
set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])

xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})

legend([hh1,hh2,hh3,hh4],'Inc. (mitigated)','HO (mitigated)','Inc. (unmitigated)','HO (unmitigated)','location','west')
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
grid on
grid minor
box on
hold off