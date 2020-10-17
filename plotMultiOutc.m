function f=plotMultiOutc(F1,X1,tvec,ddatax)
%F1 concat fb1, fb1 etc in 3rd dimension
%X1 concat x1 to x4 in second dimension
numSectors=63;
numPeriods=8;

dodiff=1;
yvar='Inc./hosp.';%'Susceptibles'; 'Hospitalisations';
solvetype=2;
tend=tvec(end);%720;%For plot only
cmap=lines(7);
cmap2=cmap+.5*(1-cmap);
col2=.5*[1,1,1];

figure
fs=10; lw=2;
%{
t2=F1(:,1,1); inc2=F1(:,2,1); h2=F1(:,3,1);
t3=F1(:,1,2); inc3=F1(:,2,2); h3=F1(:,3,2);
if dodiff==1
    inc2=-diff(inc2,1);%f1=Sout
    tdiff=diff(t2,1);
    inc2=inc2./tdiff;%repmat - older version?
    t2(1)=[];
    h2(1)=[];
    inc3=-diff(inc3,1);%f1=Sout
    tdiff=diff(t3,1);
    inc3=inc3./tdiff;%repmat - older version?
    t3(1)=[];
    h3(1)=[];
end
%}
thresh=[18057,28115,38172,53365];
titles={'H_{max}=18,057','H_{max}=28,115','H_{max}=38,172','H_{max}=53,365'};
for j=1:4
    
f1=F1(:,:,j+2);

t1=f1(:,1); inc1=f1(:,2); h1=f1(:,3);%heSimCovid19 output
if dodiff==1
    inc1=-diff(inc1,1);%f1=Sout
    tdiff=diff(t1,1);
    inc1=inc1./tdiff;%repmat - older version?
    t1(1)=[];
    h1(1)=[];
end

subplot(2,4,j)
maxY=9.5e4;%max([inc1;h1]);%inc2;h2
%
h=zeros(1,4);
hold on
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
%hh3=plot(t2,inc2,'linewidth',lw,'color',cmap2(1,:));
%hh4=plot(t2,h2,'--','linewidth',lw,'color',cmap2(2,:));
%hh5=plot(t3,inc3,'linewidth',lw,'color',cmap2(1,:));
%hh6=plot(t3,h3,'--','linewidth',lw,'color',cmap2(2,:));
hh1=plot(t1,inc1,'linewidth',lw,'color',cmap(1,:));
hh2=plot(t1,h1,'--','linewidth',lw,'color',cmap(2,:));
%}
points=tvec+10;
pointsy=.93*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',15);
text(tvec(2)+5,pointsy,'LD','fontsize',15);
for i=3:8%lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',15)
end
title(titles(j))
xlabel('Time (month)','FontSize',fs);
ylabel('Population','FontSize',fs);%yvar

xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D','J','F'})

set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])
legend([hh1,hh2],'Inc.','Hosp. occ.','location','west')
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
grid on
grid minor
box on
hold off

subplot(2,4,4+j)
if numSectors==10
    xj=reshape([ones(numSectors,1);ddatax.xmin';X1(:,j)],numSectors,numPeriods);
    colormap parula
    hold on
    %{
    bar(xj','hist')
    set(gca,'fontsize',fs)
    xticks(1:numPeriods);
    xticklabels({'PRE','LD','1','2','3','4','5','6'})
    xlabel('Intervention step (\tau)')
    ylabel('X_i(\tau)')%,'rotation',0)
    axis([.5,numPeriods+.5,0,1])
    %}
    %
    xvec=(1:numPeriods)';
    xlab=num2str(xvec-2);
    xvec=xvec-.5;
    yvec=(1:10)';
    ylab=num2str(yvec);
    yvec=yvec-.5;
    imagesc(xj)
    set(gca,'YDir','normal')
    xlabel('Period')
    ylabel('Sector')
    set(gca,'fontsize',fs,'xtick',xvec,'xticklabels',{'PRE','LD',xlab(3:end,:)},'ytick',yvec,'yticklabels',ylab)
    axis([xvec(1),xvec(end)+1,.5,10.5])%yvec(1),yvec(end)+1])
    caxis([0,1])
    colorbar
    %}
    grid on
    grid minor
    box on
elseif numSectors==63
    xj=[ones(numSectors,1);ddatax.xmin';X1(:,j)];
    %xj=X1(:,j);
    xj=reshape(xj,numSectors,numPeriods);
    xvec=(1:numPeriods)';
    xlab=num2str(xvec-2);
    xvec=xvec-.5;
    yvec=(1:10:63)';
    ylab=num2str(yvec);
    yvec=yvec-.5;
    imagesc(xj)
    set(gca,'YDir','normal')
    xlabel('Period')
    ylabel('Sector')
    set(gca,'fontsize',fs,'xtick',xvec,'xticklabels',{'PRE','LD',xlab(3:end,:)},'ytick',yvec,'yticklabels',ylab)
    axis([xvec(1),xvec(end)+1,.5,63.5])%yvec(1),yvec(end)+1])
    caxis([0,1])
    colorbar
    grid on
    box on
end
end