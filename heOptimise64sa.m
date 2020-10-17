function [xoptim,yoptim,exitflag]=heOptimise64sa%(xoptim)%(NNsector,data64,econ64)
%%
load('NNs64.mat','NNs64')
NNsector=NNs64;
load('ddata64.mat','ddata64')
load('eecon64.mat','eecon64')
%load('GS2xoutR1.mat','xoptim')
%%
hospThresh=[18057,1];%18057   28115   38172	53365

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
b2=[b+1e-6;zeros(numSect,1)];
lx=numInt*numSect;
lb=repmat(ddata64.xmin',numInt,1);%zeros(lx,1);
ub=lb;
lb(lb>1)=1;
ub(ub<1)=1;
X0=lb;
%{
%Schools open fully:
xlb=repmat(ddata64.xmin',numInt,1);
xlb((0:5)*63+55,1)=.8;
xlb(xlb>1)=1;
lb=xlb;
X0=xlb;
%X0(2*63+55:63:end)=.8;
%}
%%
fun1=@(Xit)econGDP(objFun,Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata64);
fun2=@(Xit)forSimAnn(objFun,Z2,b2,pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata64);
%
%rng default % For reproducibility
%options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'algorithm','interior-point','UseParallel',true);
%[xoptim,yoptim,exitflag,output]=simulannealbnd(fun2,X0,lb,ub);%,options);
%
%options=optimoptions('particleswarm','UseParallel',true,'UseVectorized',false);
%problem=createOptimProblem('particleswarm','x0',X0,'objective',fun2,'lb',lb,'ub',ub);
[xoptim,yoptim,exitflag,output]=particleswarm(fun2,378,X0,lb,ub);
%%
save('ps1xoutR1.mat','xoptim','yoptim','exitflag','output')
end

function f=econGDP(obj,Xit)%,lx)
f=-sum(obj.*Xit);
end

function [c,cex]=epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata64)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0,ddata64);
c=max(h-hospThresh);
cex=[];
end

function f=forSimAnn(obj,A,b,pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata64)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0,ddata64);
if max(h-hospThresh)>0
    f=0;
elseif max(A*Xit-b)>0
    f=0;
else
    f=-sum(obj.*Xit);
end
end