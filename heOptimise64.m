function [xoptim,yoptim,exitflag]=heOptimise64%(xoptim)%(NNsector,data64,econ64)
%%
load('NNs64.mat','NNs64')
NNsector=NNs64;
load('ddata64.mat','ddata64')
load('eecon64.mat','eecon64')
%%
hospThresh=[38172,1];%38172	53365	74966	82966

t0=-49.7725;%-53.8370;
t1=86.3881;%87.1771;
tend=720;
months=[1,32,61,92,122,153,183,214,245,275,306,336,361,392,420,451,481,tend];
tvec=[t0,t1,months(5:11)];%Fit tvec(1)=t0 %11
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
b2=[b+1e-6;zeros(numSect,1)];
lx=numInt*numSect;
lb=repmat(ddata64.xmin',numInt,1);%zeros(lx,1);
ub=ones(lx,1);
X0=repmat(ddata64.xmin',numInt,1);%zeros(lx,1);
%{
%Schools open fully in last 2 months:
xlb=repmat(ddata64.xmin',numInt,1);
xlb(4*63+55,1)=1;
xlb(5*63+55,1)=1;
lb=xlb;
X0=xlb;%repmat(xlb',numInt,1);%zeros(lx,1);
%}
%%
fun1=@(Xit)econGDP(objFun,Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata64);
%{
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'MaxIterations',10000000,'algorithm','interior-point');%'sqp' %'interior-point'
[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
%[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,[],options);
%}
%
rng default % For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',200000,'MaxIterations',200000,'algorithm','interior-point','UseParallel',true);%best death timeseries from public ONS
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
%problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',[],'options',options);
ms=MultiStart;
[xoptim,yoptim]=run(ms,problem,20);
exitflag='Not relevant for ms';
%save('optimOut2.mat','xoptim','yoptim','exitflag')
%}
%%
save('optimOut.mat','xoptim','yoptim','exitflag')
end

function f=econGDP(obj,Xit)%,lx)
f=-sum(obj.*Xit);
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