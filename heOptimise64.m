function [xoptim,yoptim,exitflag]=heOptimise64%(NNsector,data64,econ64)
%%
load('NNs64.mat','NNs64')
NNsector=NNs64;
load('data64.mat','data64')
load('econ64.mat','econ64')
%%
hospThresh=[10^4,1];
t0=-44.3012;
tend=720;
months=[1,32,61,92,122,153,183,214,245,275,306,336,361,392,420,451,481,tend];
tvec=[t0,months(4:11)];%Fit tvec(1)=t0 %11
lastInt=1;%(tend-tvec(end-1))/30;
numInt=length(tvec)-3;
numAges=4;
%%
%NNsector=data64.NNsector;
G=econ64.G;
b=econ64.b;
objFun=econ64.obj;%Monthly
G(:,44)=G(:,44)+G(:,45);
G(44,:)=G(44,:)+G(45,:);
G(45,:)=[];
G(:,45)=[];
b(44)=b(44)+b(45);
b(45)=[];
objFun(44)=objFun(44)+objFun(45);
objFun(45)=[];
%%
%Tune epi-midel to pre-lockdown:
[pr,NN,n,nbar,na,NNbar,NNrep,Din,beta]=hePrepCovid19(NNsector,data64);
numSect=length(NNsector)-numAges;
objFun=repmat(objFun,numInt,1);
%%
%Econ constrained over whole period:
Z2=[repmat(G,1,numInt);repmat(-G,1,numInt)];
b2=[b;zeros(numSect,1)];
lx=numInt*numSect;
lb=repmat(data64.xmin',numInt,1);%zeros(lx,1);
ub=ones(lx,1);
X0=ones(lx,1);%repmat(data64.xmin',numInt,1);%zeros(lx,1);
%%
%
fun1=@(Xit)econGDP(objFun,Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,data64);
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'algorithm','interior-point');%'sqp' %'interior-point'
[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,[],options);%nonlcon
%}
%{
rng default % For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);%best death timeseries from public ONS
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
ms=MultiStart;
[xoptim,yoptim]=run(ms,problem,20);
exitflag='Not relevant for ms';
save('optimOut2.mat','xoptim','yoptim','exitflag')
%}
%%
save('optimOut.mat','xoptim','yoptim','exitflag')
end

function f=econGDP(G,Xit)%,lx)
f=-sum(G.*Xit);
end

function [c,cex]=epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,data64)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0,data64);
c=max(h-hospThresh);
cex=[];
end
%{
%Econ constrained by month:
ball=repmat(b',numInt,1);
kkron=eye(numInt-1);
kkron(end+1,end+1)=1/lastInt;%Account for last interval being longer
Gall=sparse(kron(kkron,G));
lx=size(Gall,1);
%Gvec=diag(Gall);
Z2=[Gall;-Gall];%<0,b
b2=[ball/numInt;zeros(lx,1)];
%
lb=zeros(lx,1);
ub=ones(lx,1);
X0=.1*rand(lx,1);
%}