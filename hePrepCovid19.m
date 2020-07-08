function [pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(D,data64)%,R0)
%D - column vector or populations.
%Possible generalsiation to within-sector heterogeneity - one column per
%category. 
lc=4;
adInd=3;
lx=length(D)-lc;
%urbrur=0;%Turn in to @home vs @work
%%
%Population density:
[n,na]=size(D);
nbar=n*na;
NN=sum(D,2);
NNbar=reshape(D,nbar,1);
NNrep=repmat(NN,na,1);
%{
if urbrur==1
    hvec=kron(hvec,ones(n,1));
    muvec=kron(muvec,ones(n,1));
    %Kkron=[1,.05;.05,.75];%Sonoma
    Kkron=[1,.05;.05,.85];%Mendocino
else
    Kkron=1;
end
%}
%{
C=[.3827    0.9115    0.0419;
    1.2062    7.3538    0.1235;
    0.1459    0.4810    0.1860];
%}
%C=eye(na);

%K=heMakeDs(NN,eye(10));
D=heMakeDs(NN,ones(lx,1),data64,0);%,1);

%K=rand(n);
%K=normr(K);
%D=kron(C,K);
%%
%Hard code from heParamsAge:
Tilih=5;
Thdeath=10.81;%Median 8
Threc=12.73;%Median 9
pdeath=.39;
pili=[0.4670,0.4786,0.6590,0.7252]';
pili=[repmat(pili(adInd),lx,1);pili];
ph=[0.0090,0.0152,0.2910,0.6738]';
ph=[repmat(ph(adInd),lx,1);ph];
pd=[0.0220,0.0220,0.0416,0.3514]';
pd=[repmat(pd(adInd),lx,1);pd];

%toHosp=3;%Symp to hosp%****
Text=4.6;
Tonset=1;
pr=struct;
pr.sigma=1/Text;
pr.omega=1/Tonset;
pr.g1=1/2.1;
pr.g2=1/2.1;
pr.p1=1-.34;
pr.p2=.505*pili;%Vector
pr.p3=1;
pr.p4=1;
pr.q1=0;
pr.q2=0;
%
pr.h=ph/Tilih;%Vector
pr.gX=(1-ph)/Tilih;%Vector
pr.mu=pdeath/Thdeath;%Vector
pr.g3=(1-pd)/Threc;%Vector
%
pr.g4=1/(1/pr.g2-1/pr.q1);
if pr.q2>0
    pr.g4X=1/(1./pr.gX+1./pr.q2);%Vector
else
    pr.g4X=pr.q2*0;
end
%
pr.odds=0;
pr.qnew=0;
pr.red=2/3;
pr.R0=1.6941;%R0;%2.7106;
%%
%isdual=1 - equivalent here
%Ceff=kron(C,Ckron);%Urb/rural mixing here
Dout=D;
Deff=Dout.*repmat(NNbar,1,n*na)./repmat(NNbar',n*na,1);
%
ntot=n*na;
onesn=ones(ntot,1);
F=zeros(4*ntot,4*ntot);
%F(1:ntot,1:ntot)=Deff;
F(1:ntot,ntot+1:end)=[pr.red*Deff,repmat(Deff,1,2)];
%vvec=kron([pr.sigma;pr.g1;pr.g2;pr.gX],ones(ntot,1));
%vvec(end-ntot+1:end)=vvec(end-ntot+1:end)+pr.h;
vvec=[pr.sigma*onesn;pr.g1*onesn;pr.g2*onesn;pr.gX+pr.h];
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=diag(-pr.sigma*(1-pr.p1)*onesn);
%V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*(1-pr.p2)*onesn);
V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*(1-pr.p2));
%V(3*ntot+1:4*ntot,1:ntot)=diag(-pr.sigma*pr.p1*pr.p2*onesn);
V(3*ntot+1:4*ntot,1:ntot)=diag(-pr.sigma*pr.p1*pr.p2);
GD=F/V;
%{
%HE:
ntot=n*na;
F=zeros(7*ntot,7*ntot);
F(1:ntot,1:ntot)=Deff;%ntot+1:end)=[repmat(2/3*Deff,1,2),repmat(Deff,1,4)];
onesn=ones(ntot,1);
vvec=[pr.sigma*onesn;pr.g1*onesn;pr.omega*onesn;pr.g2*onesn;pr.g2*onesn;pr.h*onesn;pr.h*onesn];
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=diag(-pr.sigma*(1-pr.p1)*onesn);
V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*onesn);
V(3*ntot+1:4*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p3).*(1-pr.p2));
V(4*ntot+1:5*ntot,2*ntot+1:3*ntot)=diag(-pr.p3.*(1-pr.p2));
V(5*ntot+1:6*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p4).*pr.p2);
V(6*ntot+1:7*ntot,2*ntot+1:3*ntot)=diag(-pr.p4.*pr.p2);
GD=F/V;
%}
d=eigs(GD,1); R0a=max(d); beta=pr.R0/R0a;
end