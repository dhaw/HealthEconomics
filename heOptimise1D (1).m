function [xoptim,yoptim,exitflag]=heOptimise1D(thresh2)%[xoptim,yoptim,exitflag]
t0=0;
tend=720;
months=[1,32,61,92,122,153,183,214,245,275,306,336,361,392,420,451,481];
tvec=[t0,months(4:11)];%Fit tvec(1)=t0 %11
numInt=length(tvec)-3;

NNsector=[334594,3467837,2332332,7383886,1340524,1300949,373380,4031981,9741730,1828106,33704681]';
NNsector=[sum(NNsector(1:end-1));NNsector(end)];

hospThresh=[10^4,thresh2];
%%
%Tune epi-midel to pre-lockdown:
[pr,NN,n,nbar,na,NNbar,NNrep,Din,beta]=hePrepCovid19(NNsector);%,.04,.0025);
%%
lx=numInt;
lb=zeros(lx,1);
ub=ones(lx,1);
X0=rand(lx,1);
fun1=@(Xit)econGDP(Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh);
nvars=lx;
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);%'sqp' %'interior-point'
[xoptim,yoptim,exitflag]=fmincon(fun1,X0,[],[],[],[],lb,ub,nonlcon,options);
%options=optimoptions('patternsearch','InitialMeshSize',100,'UseParallel',true);
%[xoptim,yoptim,exitflag]=paretosearch(fun1,lx,[],[],[],[],lb,ub,nonlcon);%,options);
%options=optimoptions('ga','UseParallel',true);
%[xoptim,yoptim,exitflag,~]=ga(fun1,nvars,[],[],[],[],lb,ub,nonlcon,options);
%{
%rng default%For reproducibility
gs=GlobalSearch;
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'lb',lb,'ub',ub,'nonlcon',nonlcon);
[xoptim,yoptim,exitflag]=run(gs,problem);
%}
%{
rng default % For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
ms=MultiStart;
[xoptim,yoptim]=run(ms,problem,20);
exitflag=NaN;
%}
%%

end

function f=econGDP(Xit)
f=-sum(Xit);
end

function [c,cex]=epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit',tvec,0);
c=max(h-hospThresh);
cex=[];
end