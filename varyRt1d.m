function [X,Y,reasons]=varyRt1d%(NNsector,tvec)%(NNsector,X,tvec)%(c,tvec)%

%months=[1    32    61    92   122   153   183   214   245   275   306   336];
%tvec=[t0,months(4:11)];%Fit tvec(1)=t0
%NNsector=[334594,3467837,2332332,7383886,1340524,1300949,373380,4031981,9741730,1828106,33704681]';

rt=(1:.02:1.1);
numInt=6;%Number of intervention steps
lr=length(rt);
%
X=zeros(lr,numInt);
Y=zeros(lr,1);
reasons=Y;
for i=1:lr
    [xi,yi,reasoni]=heOptimise1D(rt(i));
    X(i,:)=xi;
    Y(i,:)=yi;
    reasons(i)=reasoni;
end
%}
%{
[pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19([sum(NNsector(1:10));NNsector(end)]);%If 1D
c=cell(lr,1);
for i=1:lr
    [out,~]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,X(i,:),tvec,0);
    c{i}=out;
end
%}
%{
fs=10; lw=2;
cmap=parula(lr);
cmap=flipud(cmap);
figure
hold on
%Vertical lines:
lt=length(tvec);
for i=2:lt-1
    plot(tvec(i)*[1,1],[0,10^7],'k--','linewidth',1)
end
maxY=1;
h=zeros(lr,1);
for i=1:lr
    ci=c{i};
    %if reasons(i)==2
        h(i)=plot(ci(:,1),ci(:,2),'-','linewidth',lw,'color',cmap(i,:));
        maxY=max(maxY,max(ci(:,2)));
    %else
        %h(i)=plot(ci(:,1),ci(:,2),'-','linewidth',lw,'color',.5*[1,1,1]);
    %end
end
%Numbers:
points=tvec+10;
pointsy=.93*maxY;
txt={'1','2','3','4','5','6'};
for i=3:lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',20)
end
%
set(gca,'fontsize',fs)
xlabel('Time (days)')
ylabel('Hospitalisations')
axis([tvec(1),tvec(end),0,maxY])
legend(h,num2str(rt'),'location','NW');
grid on
grid minor
box on
%}
%f=X;
%g=c;
%save('hpcRun1.mat','Xit','epiCell')