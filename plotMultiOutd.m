function f=plotMultiOutd(f1,x1,tvec,ddatax)
numPeriods=5;
numSectors=length(x1)/(numPeriods-2);
lt=length(tvec);
dodiff=1;

figure
fs=10; lw=2;
cmap=lines(2);
thresh=[12000,18000,24000];
numThresh=length(thresh);
t1=f1(:,1); inc1=f1(:,2); h1=f1(:,3);%heSimCovid19 output
if dodiff==1
    inc1=-diff(inc1,1);%f1=Sout
    tdiff=diff(t1,1);
    inc1=inc1./tdiff;%repmat - older version?
    t1(1)=[];
    h1(1)=[];
end

maxY=max([inc1;h1]);%inc2;h2
%
hold on
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for j=1:numThresh
    plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
end
hh1=plot(t1,inc1,'linewidth',lw,'color',cmap(1,:));
hh2=plot(t1,h1,'--','linewidth',lw,'color',cmap(2,:));
%}
points=tvec+10;
pointsy=.93*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',15);
text(tvec(2)+5,pointsy,'LD','fontsize',15);
for i=3:lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',15)
end
xlabel('Time','FontSize',fs);
ylabel('Number','FontSize',fs);%yvar
set(gca,'FontSize',fs);
axis ([0,tvec(end),0,maxY])

xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
if numPeriods==8
	xlabels2=({'Jan','Mar 26th','Sep','Oct','Nov','Dec','Jan','Feb'});
elseif numPeriods==5
    %xlabels2=({'Jan','Mar 26th','Sep','Nov','Jan'});
    xlabels2=({'PRE','LD','1','2','3'});
else
    error('Data missing for nunmPeriods')
end

legend([hh1,hh2],'Inc.','Hosp. occ.','location','west')
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
grid on
grid minor
box on
hold off

figure
%subplot(2,1,2)
if numSectors==10
    x1=[ones(numSectors,1);ddatax.xmin';x1];
    x=reshape(x1,numSectors,numPeriods);
    xvec=(1:numPeriods)';
    xlab=num2str(xvec-2);
    xvec=xvec-.5;
    yvec=(1:10)';
    ylab=num2str(yvec);
    yvec=yvec-.5;
    %colormap gray
    imagesc([0,tvec(2:end)],0:numSectors,x)
    set(gca,'YDir','normal')
    xlabel('Period')
    ylabel('Sector')
    set(gca,'fontsize',fs,'xtick',[0,tvec(2:end-1)],'xticklabels',xlabels2,'ytick',yvec,'yticklabels',ylab);%{'PRE','LD',xlab(3:end,:)}
    axis([xvec(1),xvec(end)+1,.5,10.5])%yvec(1),yvec(end)+1])
    caxis([0,1])
    colorbar
    grid on
    box on
elseif numSectors==63
    x1=[ones(numSectors,1);ddatax.xmin';x1];
    x=reshape(x1,numSectors,numPeriods);
    xvec=(1:numPeriods)';
    xlab=num2str(xvec-2);
    xvec=xvec-.5;
    yvec=(1:10:63)';
    ylab=num2str(yvec);
    yvec=yvec-.5;
    %colormap gray
    imagesc(x)
    set(gca,'YDir','normal')
    xlabel('Period')
    ylabel('Sector')
    set(gca,'fontsize',fs,'xtick',[0,tvec(2:end-1)],'xticklabels',xlabels2,'ytick',yvec,'yticklabels',ylab);%{'PRE','LD',xlab(3:end,:)}
    axis tight%([0,numPeriods,.5,63.5])%yvec(1),yvec(end)+1])
    caxis([0,1])
    colorbar
    grid on
    box on
end