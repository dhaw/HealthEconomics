function f=heOptimise1D
t0=0;
tend=720;
months=[1    32    61    92   122   153   183   214   245   275   306   336];
tvec=[t0,months(4:11)];%Fit tvec(1)=t0
numInt=length(tvec)-3;

NNsector=[334594,3467837,2332332,7383886,1340524,1300949,373380,4031981,9741730,1828106,33704681]';
NNsector=[sum(NNsector(1:end-1));NNsector(end)];

hospThresh=[.5*10^5,1];
%%
%Tune epi-midel to pre-lockdown:
[pr,NN,n,nbar,na,NNbar,NNrep,Din,beta]=hePrepCovid19(NNsector);%,.04,.0025);
%%
lx=numInt;
lb=zeros(lx,1);
ub=ones(lx,1);
X0=ones(lx,1);

fun1=@(Xit)econGDP(Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh);
options=optimoptions(@fmincon,'UseParallel',1,'MaxFunctionEvaluations',100000,'algorithm','interior-point');
xoptim=fmincon(fun1,X0,[],[],[],[],lb,ub,nonlcon,options);
f=xoptim;
%%

end

function f=econGDP(Xit)
f=-sum(Xit);
end

function [c,cex]=epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0);
c=max(h-hospThresh);
cex=[];
end