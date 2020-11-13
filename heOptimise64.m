function [xoptim,yoptim,exitflag]=heOptimise64%(xoptim)%(NNsector,data64,econ64)
%%
%Scenario (A or B):
scenB=1;

%Intervention intervals (1 for 6x1, 2 for 3x2):
monthPeriod=2;

%Output file name:
fileName='gs_B2.mat';

%Option to input X0:
X0in=0;
%Note: if numPeriod=1, code assumes that X0 is a 3x2 solution

%Epi constraints (H_max, R_end):
hospThresh=[18000,1];
%%
load('NNs64.mat','NNs64')
NNsector=NNs64;
load('ddata64.mat','ddata64')
load('eecon64.mat','eecon64')
if X0in==1
    %Input XO as 'xoptim in this file:
    load('gs_1.mat','xoptim')
end
%%
t0=-49.7725;%-49.7725;
t1=86.3881;%86.3881;
tend=720;
months=[1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,tend];
tvec=[t0,t1,months(9:monthPeriod:15)];%Fit tvec(1)=t0 %11
lastInt=1;%(tend-tvec(end-1))/30;
numInt=length(tvec)-3;
numAges=4;
%%
%NNsector=data64.NNsector;
G=eecon64.G*monthPeriod;
b=eecon64.b;
objFun=eecon64.obj*monthPeriod;%Monthly
%%
%Tune epi-midel to pre-lockdown:
[pr,NN,n,nbar,na,NNbar,NNrep,Din,beta]=hePrepCovid19(NNsector,ddata64);
numSect=length(NNsector)-numAges;
objFun=repmat(objFun,numInt,1);
%%
%Econ constrained over whole period:
Z2=[repmat(G,1,numInt);repmat(-G,1,numInt)];
%b2=[b+1e-6;zeros(numSect,1)];
b2=[b;zeros(numSect,1)]+4000;

Z2=repmat(-G,1,numInt);
b2=max(-6*G*ddata64.xmin',0);

lx=numInt*numSect;
lb=repmat(ddata64.xmin',numInt,1);%zeros(lx,1);
ub=lb;
lb(lb>1)=1;
ub(ub<1)=1;
%%
if X0in==0
    X0=lb;
else
    X0=xoptim;
    if monthPeriod==1
        X0=[X0(1:63);X0(1:63);X0(64:126);X0(64:126);X0(127:189);X0(127:189)];
    end
end
%Chunk at bottom
fun1=@(Xit)econGDP(objFun,Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata64);
if X0in==0
    if scenB==1
        X0(55:63:end)=0.8;
    end
    X0=linprog(ones(1,size(Z2,2)),Z2,b2,[],[],X0,X0+.1*(ub-X0));
end
%lb=zeros(378,1);

%Set lower bound:
lb=0.8*lb;

%Health:
%lb(56:63:end)=ddata64.xmin(56);
if scenB==1
	lb(55:63:end)=0.8;
end
%%
%{
options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'algorithm','interior-point');%'sqp' %'interior-point'
[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
%[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,[],options);
%}
%%
%{
rng default % For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
ms=MultiStart('StartPointsToRun','bounds');%,'Display','iter','MaxTime',3600);
%rs=RandomStartPointSet('NumStartPoints',3);
%allpts={rs};
%[x fval eflag output manymins]=run(ms,problem,allpts);
%
%numPoints=12;
%sp=repmat(lb',numPoints,1)+repmat(ub'-lb',numPoints,1).*rand(numPoints,lx);
%sp=CustomStartPointSet(sp);
[xoptim,yoptim,exitflag,output,manymins]=run(ms,problem);%,sp);
%exitflag='Not relevant for ms';
%}
%%
%{
rng default
options=optimoptions('patternsearch','UseParallel',true,'UseVectorized', false,'UseCompletePoll',true);%,'InitialPoints',ips);%'MaxFunctionEvaluations',1e6,'MaxIterations',1e6);%,'ConstraintTolerance',1);
[xoptim,yoptim,exitflag]=patternsearch(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
%}
%%
%{
rng default
%ips=struct;
%ips.X0=[.9*X0';.8*X0'+.5*(ub'-X0')];
np=23;
xmin=repmat(.9*X0',np,1);
xmax=repmat(.1*X0'+.5*(ub'-X0'),np,1);
ips=xmin+rand(np,378/monthPeriod).*xmax;
ips=[X0';ips];
fun2=@(Xit)econGDPga(objFun,Xit);
optimoptions('patternsearch','UseVectorized',true);%,'UseCompletePoll',true);
options=optimoptions('paretosearch','InitialPoints',ips);%'UseParallel',false,'UseVectorized', true,'InitialPoints',ips);%,'ParetoSetSize',1e4 ,'UseCompletePoll',true
[xoptim,yoptim,exitflag]=paretosearch(fun2,378/monthPeriod,Z2,b2,[],[],lb,ub,nonlcon,options);
%}
%%
%{
fun2=@(Xit)econGDPga(objFun,Xit);
rng default
options=optimoptions('ga','UseParallel',true,'UseVectorized',false,'InitialPopulationRange',[X0';X0'+.1*(ub'-X0')]);
%opts.InitialPopulationRange=[lb'; X0'];%2 rows, nvar columns - lb;ub
[xoptim,yoptim,exitflag]=ga(fun2,378/monthPeriod,Z2,b2,[],[],lb,ub,nonlcon,options);
xoptim=xoptim';
%}
%%
%
rng default % For reproducibility
options=optimoptions(@fmincon,'UseParallel',true,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6);%,'algorithm','interior-point',
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
gs=GlobalSearch;
[xoptim,yoptim,exitflag]=run(gs,problem);
%}
%%
save(fileName,'xoptim','yoptim','exitflag')%,'output','manymins')
end

function f=econGDP(obj,Xit)%,lx)
f=-sum(obj.*Xit);
end

function f=econGDPga(obj,Xit)%,lx)
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
%%
%X0=repmat(ddata64.xmin',numInt,1);%zeros(lx,1);
%X0=[ones(2*63,1);repmat(ddata64.xmin',4,1)];
%xunder=ddata64.xmin'+.4*(1-ddata64.xmin');
%X0=[repmat(xunder,2,1);repmat(ddata64.xmin',4,1)];
%X0(end-62:end)=ddata64.xmin';
%
%{
%Schools open fully/to given degree:
xlb=repmat(ddata64.xmin',numInt,1);
xlb((0:5)*63+55,1)=1;
xlb(xlb>1)=1;
lb=xlb;
X0=xlb;
%}
%X0(X0<.7)=.7;
%X0(end-62:end)=ddata64.xmin';
%%
%X0(2*63+55:63:end)=.8;
%Start-point based on existing output:
%X0=xoptim;
%xmin64=repmat(ddata64.xmin',6,1);
%X0=xmin64+.8*(xoptim-xmin64);
%X0(5*63+1:end)=xmin64(5*63+1:end);
%X0=xmin64+kron([.5;.4;.3;.2;.1;0],ones(63,1)).*(xoptim-xmin64);
%Schox0(x0<.7)=.7;ols open:
%X0(4*63+55,1)=1; X0(5*63+55,1)=1;