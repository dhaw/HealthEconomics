function [f,g]=heSimCovid19vax(pr,beta,tvec,Dvec,n,nbar,NNvec,phi1,phi2,seedvec,S0,tau,plotTau)
hospInc=0;
lx=length(S0)-4;
adInd=3;
%%
%Vacconation parameters:
vx=struct;
NNage=[4064198,12192593,36577778,13005432]';%NumAgeGroupsExplicit
tpoints=[336,367,457];%Period end points - 1st Dec, 1st Jan, 1st April
tint=diff(tpoints);%Period lengths
%PERIOD 1
vx.startv1=92;%336;%336 - 1st Dec
vx.rate1=[zeros(nbar-1,1);2e6/tint(1)]/NNage(4);
vx.uptake=ones(nbar,1);
%vx.uptake1=ones(nbar,1);
%PERIOD 2:
%Proportion of working age over 50:
vx.startv2=102;%367;%Next vax change; 398=1st Feb 2021
prop50=15/45;
%vx.uptake2=;
vx.startv2=102;%367;%Next vax change; 398=1st Feb 2021
vx.rate2=[prop50*ones(nbar-4,1);0;0;prop50;1]/tint(2);%*NNnext/NNnext
%PERIOD 3:
%
%Fix parameters:
%VACCINE 1:
vx.sus1=.9;
vx.tr1=0;
vx.p1v1=2/3;
vx.p2v1=pr.p2;
%VACCINE 2:
vx.sus2=.9;
vx.tr2=0;
vx.p1v2=2/3;
vx.p2v2=pr.p2;
%%
solvetype=2;
zn=zeros(nbar,1);
NNbar=NNvec(:,1);
if solvetype==2
    lt=length(tvec);
    t0=tvec(1);%tvec(1)=start time for while simulation
    %y0=[S0;repmat(zn,11,1);NNbar-S0];%HE
    y0=[S0;repmat(zn,18,1);NNbar-S0];%IC
    toutAll=[];
    Sout=[];
    Soutv1=Sout;
    Soutv2=Sout;
    Hout=Sout;
    Dout=Sout;
    Iout=[];
    DEout=zeros(nbar,lt);
    Rout=DEout;
    Rt=zeros(lt-1,1);
    HnewAll=[];
    for i=1:lt-1
        D=Dvec(:,:,i);
        tend=tvec(i+1);
        NNfeed=NNvec(:,i);
        
        %vx.rate1=2e6*NNfeed/sum(NNfeed).*ones(nbar,1);
        
        NNfeed(NNfeed==0)=1;
        %
        %pr.z=0;
        %{
        if i>=2%3
            pr.q1=pr.qnew;%Won't turn off****
            pr.q2=pr.qnew;
            %pr.z=pr.znew;
            pr.odds=pr.oddsnew;
        end
        %}
        %topen=1;%tvec(3);
        if hospInc>=1
            Rt(i)=heComputeEigs(pr,beta,D,NNfeed,nbar,y0(1:lx+4));
        end
        [tout,Sclass,Sclassv1,Sclassv2,Hclass,Dclass,DEcum,Rcum,Itot,y0,Hnew]=integr8(pr,beta,nbar,NNfeed,D,phi1,phi2,seedvec,t0,tend,y0,i,hospInc,vx);
        toutAll=[toutAll;tout(2:end)];
        Sout=[Sout;Sclass(2:end,:)];
        Soutv1=[Soutv1;Sclassv1(2:end,:)];
        Soutv2=[Soutv2;Sclassv2(2:end,:)];
        Hout=[Hout;Hclass(2:end,:)];
        Dout=[Dout;Dclass(2:end,:)];
        DEout(:,i)=DEcum;
        HnewAll=[HnewAll;Hnew(2:end,:)];
        %Rout(:,i)=Rcum;
        Iout=[Iout;Itot(2:end,:)];
        t0=tend;
        %%
        %HE:
        %Calculate in advance?
        %Need to add age structure in here - only uses n, not na****
        %
        %Rt(i)=heComputeEigs(pr,beta,D,NNfeed,nbar,Sclass(end,:)');%****
        if hospInc==0
            %Rt(i)=heComputeEigs(pr,beta,D,NNfeed,nbar,y0(1:lx+4));
            Rt(i)=heComputeEigs(pr,beta,D,NNfeed,nbar,Sclass(end,:)');
        end
        %
        if i<lt-1
            Xh2w=NNvec(1:lx,i+1)-NNvec(1:lx,i);%Addition to each wp next intervention step
            Xw2h=-Xh2w; Xw2h(Xw2h<0)=0; 
            Xw2h=Xw2h./NNvec(1:lx,i);
            Xw2h(NNvec(1:lx,i)==0)=0;
            if NNvec(lx+adInd,i)>0
                Xh2w(Xh2w<0)=0;
                Xh2w=Xh2w/NNvec(lx+adInd,i);
            else
                Xh2w=0;
            end
            %Move all infection statuses:
            %y0=reshape(y0,[n,13]);%HE
            y0=reshape(y0,[n,20]);%IC
            %y0w2h=y0(1:n-1,:).*repmat(Xw2h,1,13);%HE
            y0w2h=y0(1:lx,:).*repmat(Xw2h,1,20);%IC
            y0w2h=[-y0w2h;sum(y0w2h,1)];
            y0h2w=y0(lx+adInd,:);
            y0h2w=kron(y0h2w,Xh2w);
            y0h2w=[y0h2w;-sum(y0h2w,1)];
            y0([1:lx,lx+adInd],:)=y0([1:lx,lx+adInd],:)+y0w2h+y0h2w;
            %y0=reshape(y0,13*nbar,1);%HE
            y0=reshape(y0,20*nbar,1);%IC
        end
        %%
        %Calculate per-person vaccination rates
    end
    if tau==plotTau
            dodiff=1;
            plotEpi(toutAll,Sout+Soutv1+Soutv2,Hout,n,dodiff,tvec);%-tvec(2)+38 %Sout+Soutv1+Soutv2
    end
elseif solvetype==1
    error('Final size calculations not possible')
elseif solvetype==3
    error('Code not written yet')
    %f=stochSim(y0,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau,alpha);
end
%%
%Function outputs here:
%For plots:

if hospInc==0
    %Rt calculated at end of each period
    %f=[toutAll(toutAll>0),sum(Hout(toutAll>0,:),2)];%For epi fit
    f=[toutAll(toutAll>0),sum(Sout(toutAll>0,:),2),sum(Hout(toutAll>0,:),2)];%sum(DEout(end,:));
    g=sum(sum(Hout(toutAll>tvec(3),:)));%[max(sum(Hout(toutAll>tvec(3),:),2)),Rt(end)];%Main constraints %DEout;
elseif hospInc==1%****
    f=[toutAll(toutAll>0),HnewAll(toutAll>0,:)];
    g=sum(sum(HnewAll(toutAll>tvec(3),:)));
    %g=Rt(3);
elseif hospInc==2
    f=[toutAll(toutAll>0),sum(Sout(toutAll>0,:),2),sum(Hout(toutAll>0,:),2)];%Occupancy
    g=sum(sum(HnewAll(toutAll>tvec(3),:)));%Incidence
end
end

function [tout,Sclass,Sclassv1,Sclassv2,Hclass,Dclass,DEcum,Rcum,Itot,y0new,Hnew]=integr8(pr,beta,nbar,NN0,D,phi1,phi2,seedvec,t0,tend,y0,topen,hospInc,vx)
%ncomps=13;%Number of compartments

    %varPassedOut=0;
    
    fun=@(t,y)integr8covid(t,y,pr,beta,nbar,NN0,D,phi1,phi2,seedvec,topen,vx);
    [tout,yout]=ode45(fun,(round(t0):1:tend),y0);
    %%
    %
    %For fit:
    if hospInc==1
        Hnew=zeros(length(tout),nbar);
        for i=1:length(tout)
            [~,Hnewi]=fun(tout(i),yout(i,:)');
            Hnew(i,:)=Hnewi';
        end
    else
        Hnew=0;
    end
    %%
    Sclass=yout(:,1:nbar);
    Sclassv1=yout(:,5*nbar+1:6*nbar)+yout(:,6*nbar+1:7*nbar);
    Sclassv2=yout(:,11*nbar+1:12*nbar)+yout(:,12*nbar+1:13*nbar);
    %{
    %HE:
    Hclass=yout(:,10*nbar+1:11*nbar);
    Dclass=yout(:,11*nbar+1:12*nbar);
    DEcum=yout(end,11*nbar+1:12*nbar);%Deaths
    Rcum=yout(end,12*nbar+1:end);
    Itot=sum(yout(:,4*nbar+1:8*nbar),2);
    %}
    %
    %IC:
    Hclass=yout(:,17*nbar+1:18*nbar);
    Dclass=yout(:,18*nbar+1:19*nbar);
    DEcum=yout(end,18*nbar+1:19*nbar);%Deaths
    Rcum=yout(end,19*nbar+1:end);
    Itot=sum(yout(:,2*nbar+1:5*nbar),2)+sum(yout(:,8*nbar+1:11*nbar),2)+sum(yout(:,14*nbar+1:17*nbar),2);%Wrong in non-vax code but unused
    y0new=yout(end,:)';
    %}
end

function [f,g]=integr8covid(t,y,pr,betaIn,nbar,NN0,D,phi1,phi2,seedvec,itime,vx)
%phi=phi1-phi2*cos(pi*t/180);%Seasonality****
phi=phi1;
%%
S=y(1:nbar);
E=y(nbar+1:2*nbar);
Ia=y(2*nbar+1:3*nbar);
Im=y(3*nbar+1:4*nbar);
Is=y(4*nbar+1:5*nbar);
%
Sholdv1=y(5*nbar+1:6*nbar);
Sv1=y(6*nbar+1:7*nbar);
Ev1=y(7*nbar+1:8*nbar);
Iav1=y(8*nbar+1:9*nbar);
Imv1=y(9*nbar+1:10*nbar);
Isv1=y(10*nbar+1:11*nbar);
%
Sholdv2=y(11*nbar+1:12*nbar);
Sv2=y(12*nbar+1:13*nbar);
Ev2=y(13*nbar+1:14*nbar);
Iav2=y(14*nbar+1:15*nbar);
Imv2=y(15*nbar+1:16*nbar);
Isv2=y(16*nbar+1:17*nbar);
%
H=y(17*nbar+1:18*nbar);
DE=y(18*nbar+1:19*nbar);
R=y(19*nbar+1:20*nbar);
I=pr.red*Ia+Im+Is+vx.tr1*(pr.red*Iav1+Imv1+Isv1)+vx.tr2*(pr.red*Iav2+Imv2+Isv2);
seed1=seedvec.*S./NN0;%No seed in sus vaccinated
beta=betaIn;
%{
if t>topen && t<topen+30
    beta=betaIn*(.25+.75*(t-topen)/30);
else
    beta=betaIn;
end
%}
foi=phi*beta*(D*(I./NN0));
Sfoiv1=phi*((1-vx.sus1)*beta*Sv1.*(D*(I./NN0)));%+seed1);
Sfoiv2=phi*((1-vx.sus2)*beta*Sv2.*(D*(I./NN0)));%+seed1);
%
%Uptake
if t>vx.startv2
    vrate1=.2*vx.rate2.*S;
    vrate2=.8*vx.rate2.*S;
elseif t>vx.startv1
    vrate1=vx.rate1.*S;
    vrate2=zeros(nbar,1);
else 
    vrate1=zeros(nbar,1);
    vrate2=zeros(nbar,1);
end
%
Sdot=-S.*foi-phi*seed1-vrate1-vrate2;
Edot=(S+Sholdv1+Sholdv2).*foi+phi*seed1-pr.sigma*E;
Iadot=(1-pr.p1)*pr.sigma*E-pr.g1*Ia;
Imdot=pr.p1*pr.sigma*(1-pr.p2)*pr.p3.*E-pr.g2.*Im;
Isdot=pr.p1*pr.sigma*pr.p2*pr.p4.*E-(pr.h+pr.gX).*Is;
%
Sholddotv1=vrate1-Sholdv1.*foi-Sholdv1/28;
Sdotv1=Sholdv1/28-Sfoiv1;
Edotv1=Sfoiv1-pr.sigma*Ev1;
Iadotv1=(1-vx.p1v1)*pr.sigma*Ev1-pr.g1*Iav1;
Imdotv1=vx.p1v1*pr.sigma*(1-vx.p2v1)*pr.p3.*Ev1-pr.g2.*Imv1;
Isdotv1=vx.p1v1*pr.sigma*vx.p2v1*pr.p4.*Ev1-(pr.h+pr.gX).*Isv1;
%
Sholddotv2=vrate2-Sholdv2.*foi-Sholdv2/28;
Sdotv2=Sholdv2/28-Sfoiv2;
Edotv2=Sfoiv2-pr.sigma*Ev2;
Iadotv2=(1-vx.p1v2)*pr.sigma*Ev2-pr.g1*Iav2;
Imdotv2=vx.p1v2*pr.sigma*(1-vx.p2v2)*pr.p3.*Ev2-pr.g2.*Imv2;
Isdotv2=vx.p1v2*pr.sigma*vx.p2v2*pr.p4.*Ev2-(pr.h+pr.gX).*Isv2;
%
Hdot=pr.h.*(Is+Isv1+Isv2)-(pr.g3+pr.mu).*H;
DEdot=pr.mu.*H;
Rdot=pr.g1*(Ia+Iav1+Iav2)+pr.g2.*(Im+Imv1+Imv2)+pr.g3.*H+pr.gX.*(Is+Isv1+Isv2);
f=[Sdot;Edot;Iadot;Imdot;Isdot;Sholddotv1;Sdotv1;Edotv1;Iadotv1;Imdotv1;Isdotv1;Sholddotv2;Sdotv2;Edotv2;Iadotv2;Imdotv2;Isdotv2;Hdot;DEdot;Rdot];
g=pr.h.*(Is+Isv1+Isv2);%Hin
%Track vaccine admin to date
end

function f=plotEpi(tout,Y,H,n,dodiff,tvec)
yvar='Inc./hosp.';%'Susceptibles'; 'Hospitalisations';
solvetype=2;
tend=tvec(end);%720;%For plot only
na=size(Y,2)/n;
cmap=lines(7);
if solvetype==2
    if dodiff==1
        Y=-diff(Y,1);
        tdiff=diff(tout,1);
        Y=Y./repmat(tdiff,1,n*na);%repmat - older version?
        tout(1)=[];
        H(1,:)=[];
    end
    figure
    fs=10; lw=2;
    if na==4
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2),sum(Y(:,3*n+1:end),2)];
        Hall=[sum(H(:,1:n),2),sum(H(:,n+1:2*n),2),sum(H(:,2*n+1:3*n),2),sum(H(:,3*n+1:end),2)];
        %Yall=Y(:,1:n)+Y(:,n+1:2*n)+Y(:,2*n+1:3*n)+Y(:,3*n+1:end);
    elseif na==5
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2),sum(Y(:,3*n+1:4*n),2),sum(Y(:,4*n+1:end),2)];
    elseif na==3
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2)];
        Hall=[sum(H(:,1:n),2),sum(H(:,n+1:2*n),2),sum(H(:,2*n+1:3*n),2)];
    elseif na==1
        Yall=sum(Y,2);
        Hall=sum(H,2);
    else
        error('Number of age groups not recognised for plotting')
    end
    %maxY=max(max(Yall));
    maxY=max(max(max(Yall)),max(max(Hall)));
    %
    %Unlogged plots:
    h=zeros(1,na);
    hold on
    for i=2:length(tvec)-1
        plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
    end
    for i=1:na
        %h(i)=plot(tout,Yall(:,i),'linewidth',lw,'color',cmap(i,:));
        h1=plot(tout,Yall(:,i),'linewidth',lw,'color',cmap(1,:));%(i,:)
        h2=plot(tout,Hall(:,i),'-','linewidth',lw,'color',cmap(2,:));%'--', (i,:)
    end
    %plot(tout,Yall,'linewidth',.5,'color',cmap(2,:));
    %plot(tout,sum(Yall,2),'linewidth',lw,'color','k');
    %}
    %{
    %Logged plots:
    hold on
    semilogy(tout,Yall);
    %}
    lt=length(tvec);
    points=[0,tvec(2:end)]+10;
    pointsy=.93*maxY;
    txt={'PRE','LD','1','2','3','4','5','6'};
    for i=1:lt-1
        text(points(i),pointsy,txt{i},'fontsize',20)
    end
    
    xlabel('Time','FontSize',fs);
    ylabel('Population','FontSize',fs);%yvar
    set(gca,'FontSize',fs);
    axis([0,tend,0,maxY])
    
    xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
    xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
    
    legend([h1,h2],'Inc.','HO','location','W')
    %
    if na==4
        legend(h,{'0-4','5-19','20-64','65+'},'location','NE')
    elseif na==5
        legend(h,{'0-4','5-17','18-49','50-64','65+'},'location','NE')
    elseif na==3
        legend(h,{'0-15','16-64','65+'},'location','NE')
    end
    %}
    grid on
    grid minor
    box on
    hold off
%}
end
end