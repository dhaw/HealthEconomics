function [xoptim,yoptim,exitflag]=heOptimise64hold%(xoptim)%(NNsector,data64,econ64)
%%
load('NNs64.mat','NNs64')
NNsector=NNs64;
load('ddata64.mat','ddata64')
load('eecon64.mat','eecon64')
%load('GS2xoutR1.mat','xoptim')
%%
hospThresh=[18057,1];%38172	53365	74966	82966 % 18057   28115   38172	53365

t0=-49.7725;%-53.8370;
t1=86.3881;%87.1771;
tend=720;
months=[1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,tend];
tvec=[t0,t1,months(9:15)];%Fit tvec(1)=t0 %11
lastInt=1;%(tend-tvec(end-1))/30;
numInt=length(tvec)-3;
numAges=4;
%%
%NNsector=data64.NNsector;
G=eecon64.G;
b=eecon64.b;
objFun=eecon64.obj;%Monthly
%%
%Tune epi-midel to pre-lockdown:
[pr,NN,n,nbar,na,NNbar,NNrep,Din,beta]=hePrepCovid19(NNsector,ddata64);
numSect=length(NNsector)-numAges;
objFun=repmat(objFun,numInt,1);
%%
%Econ constrained over whole period:
Z2=[repmat(G,1,numInt);repmat(-G,1,numInt)];
%b2=[b+1e-6;zeros(numSect,1)];
b2=[b;zeros(numSect,1)];
lx=numInt*numSect;
lb=repmat(ddata64.xmin',numInt,1);%zeros(lx,1);
ub=lb;
lb(lb>1)=1;
ub(ub<1)=1;
%
X0=lb;
%X0=repmat(ddata64.xmin',numInt,1);%zeros(lx,1);
%X0=[ones(2*63,1);repmat(ddata64.xmin',4,1)];
%xunder=ddata64.xmin'+.4*(1-ddata64.xmin');
%X0=[repmat(xunder,2,1);repmat(ddata64.xmin',4,1)];
%X0=xoptim;
%X0(end-62:end)=ddata64.xmin';
%
%{
%Schools open fully:
xlb=repmat(ddata64.xmin',numInt,1);
xlb((0:5)*63+55,1)=.8;
xlb(xlb>1)=1;
lb=xlb;
X0=xlb;
%X0(2*63+55:63:end)=.8;
%}
%Start-point based on existing output:
%X0=xoptim;
%xmin64=repmat(ddata64.xmin',6,1);
%X0=xmin64+.8*(xoptim-xmin64);
%X0(5*63+1:end)=xmin64(5*63+1:end);
%X0=xmin64+kron([.5;.4;.3;.2;.1;0],ones(63,1)).*(xoptim-xmin64);
%Schools open:
%X0(4*63+55,1)=1; X0(5*63+55,1)=1;
%%
fun1=@(Xit)econGDP(objFun,Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata64);

X0=linprog(ones(1,size(Z2,2)),Z2,b2,[],[],lb,lb+.1*(ub-lb));

%%
%{
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'MaxIterations',10000000,'algorithm','interior-point');%'sqp' %'interior-point'
[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
%[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,[],options);
%}
%%
%{
rng default % For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
%problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',[],'options',options);
%ms=MultiStart;
ms=MultiStart('StartPointsToRun','bounds');%,'Display','iter','MaxTime',3600);
%rs=RandomStartPointSet('NumStartPoints',3);
%allpts={rs};
%[x fval eflag output manymins]=run(ms,problem,allpts);
numPoints=12;
sp=repmat(X0',numPoints,1)+repmat(ub'-X0',numPoints,1).*rand(numPoints,lx);
sp=CustomStartPointSet(sp);
[xoptim,yoptim,exitflag,output,manymins]=run(ms,problem,sp);
%exitflag='Not relevant for ms';
%}
%%
%{
rng default
ips=struct;
ips.X0=[X0';X0'+.1*(ub'-X0')];
options=optimoptions('patternsearch','UseParallel',true,'UseVectorized', false,'UseCompletePoll',true,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6,'ConstraintTolerance',1);
[xoptim,yoptim,exitflag]=patternsearch(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
%%options=optimoptions('paretosearch','UseParallel',true,'UseVectorized', false,'InitialPoints',ips);%,'UseCompletePoll',true
[xoptim,yoptim,exitflag]=paretosearch(fun1,378,Z2,b2,[],[],lb,ub,nonlcon,options);
%}
%%
%
rng default
options=optimoptions('ga','UseParallel',true,'UseVectorized',false,'InitialPopulationRange',[X0';X0'+.1*(ub'-X0')]);
%opts.InitialPopulationRange=[lb'; X0'];%2 rows, nvar columns - lb;ub
[xoptim,yoptim,exitflag]=ga(fun1,378,Z2,b2,[],[],lb,ub,nonlcon,options);
%}
%%
%{
rng default % For reproducibility
options=optimoptions(@fmincon,'UseParallel',true,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6);%,'algorithm','interior-point',
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
gs=GlobalSearch;
[xoptim,yoptim,exitflag]=run(gs,problem);
%}
%%
save('ga2xoutR1.mat','xoptim','yoptim','exitflag')%,'output','manymins')
end

function f=econGDP(obj,Xit)%,lx)
f=-sum(obj'.*Xit);
end

function [c,cex]=epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata64)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0,ddata64);
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
%{
%G(:,44)=G(:,44)+G(:,45);
%G(44,:)=G(44,:)+G(45,:);
G(45,:)=[];
G(:,45)=[];
%b(44)=b(44)+b(45);
b(45)=[];
%objFun(44)=objFun(44)+objFun(45);
objFun(45)=[];
%
%data64.xmin=zeros(1,63);
b(2)=0;%b(1);
%
objFun(2)=objFun(1);
g1r=G(1,:);
gc1=G(:,1);
G(2,:)=g1r;
G(:,2)=gc1;
%}