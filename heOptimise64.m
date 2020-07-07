function [xoptim,yoptim,exitflag]=heOptimise64(thresh2)
hospThresh=[10^4,thresh2];
t0=0;
tend=720;
months=[1,32,61,92,122,153,183,214,245,275,306,336,361,392,420,451,481];
tvec=[t0,months(4:11)];%Fit tvec(1)=t0 %11
lastInt=1;%(tend-tvec(end-1))/30;
numInt=length(tvec)-3;
numAges=4;
%%
NNsector=[];
G=[];
b=[];
objFun=[];%Monthly
%%
%Tune epi-midel to pre-lockdown:
[pr,NN,n,nbar,na,NNbar,NNrep,Din,beta]=hePrepCovid19(NNsector);
numSect=length(NNsector)-numAges;
objFun=repmat(objFun',numInt,1);
%%
%Econ constrained over whole period:
Z2=[repmat(G,1,numInt);repmat(-G,1,numInt)];
b2=[b';zeros(numSect,1)];
lx=numInt*numSect;
lb=zeros(lx,1);
ub=ones(lx,1);
X0=zeros(lx,1);
%%
%
fun1=@(Xit)econGDP(objFun,Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh);
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'algorithm','interior-point');%'sqp' %'interior-point'
[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
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

end

function f=econGDP(G,Xit)%,lx)
f=-sum(G.*Xit);
end

function [c,cex]=epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0);
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