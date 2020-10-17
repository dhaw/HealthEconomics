function f=heOptimise1D
load('NNs1.mat','NNs1')
NNsector=NNs1;
load('ddata1.mat','ddata1')
%
hospThresh=[38172,1];%38172	53365   74966	82966
t0=-49.7725;%-44.9779;
t1=86.3881;%86.3881;
tend=720;
months=[1,32,61,92,122,153,183,214,245,275,306,336,361,392,420,451,481,tend];
tvec=[t0,t1,months(5:11)];%Fit tvec(1)=t0 %11
%lastInt=1;%(tend-tvec(end-1))/30;
numInt=length(tvec)-3;
%numAges=4;
%{
%NNsector=[334594,3467837,2332332,7383886,1340524,1300949,373380,4031981,9741730,1828106,33704681]';
%NNsector=[sum(NNsector(1:end-1));NNsector(end)];
lx=10; lc=3;
NNsector=[334594,3467837,2332332,7383886,1340524,1300949,373380,4031981,9741730,1828106,4064198,12192593,4442459,13005432]';
NNsector=[sum(NNsector(1:lx));NNsector(lx+1:lx+lc)];
%}
%%
%Tune epi-model to pre-lockdown:
[pr,NN,n,nbar,na,NNbar,NNrep,Din,beta]=hePrepCovid19(NNsector,ddata1);%,.04,.0025);
%%
lb=zeros(numInt,1);%repmat(data1.xmin,numInt,1);
ub=ones(numInt,1);
X0=repmat(ddata1.xmin,numInt,1);%zeros(lx,1);
%Xit=repmat(ddata1.xmin,numInt,1);
%%
fun1=@(Xit)econGDP(Xit);
%fun2=@(pmod)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,data1,pmod);
nonlcon=@(pmod)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,X0,tvec,hospThresh,ddata1);
%{
options=optimoptions(@fmincon,'UseParallel',false,'MaxFunctionEvaluations',100000,'algorithm','interior-point');
xoptim=fmincon(fun1,X0,[],[],[],[],lb,ub,nonlcon,options);
%}
rng default % For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);%best death timeseries from public ONS
%problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'lb',lb,'ub',ub,'nonlcon',nonlcon2,'options',options);
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
ms=MultiStart;
[xoptim,~]=run(ms,problem,20);
%exitflag='Not relevant for ms';
save('optimOut1D.mat','xoptim')
%%
f=xoptim;
end

function f=econGDP(Xit)
f=-sum(Xit);
end

function [c,cex]=epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,data1)%,pmod)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0,data1);%,pmod);
c=max(h-hospThresh);
cex=[];
end