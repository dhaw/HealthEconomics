function f=fitEpi(NNsectorAges,ydata,data64)
%Data=matrix - vector of weekly deaths by week of year
%hePrepCovid19 - feed in R0
%heSimCovid19 - output deaths
wend=11;
xdata=(1:wend)';
%ydata=cumsum(ydata(1:wend));
v7=(5:7:wend*7-2)';
%{
x0=[2.5,0];%R0,t0
lb=[1,-60];
ub=[4,60];
%}
x0=0;
lb=0;%-100;
ub=1;%100;

fun=@(params,xdata)sim2fit(NNsectorAges,params,v7,data64);
%poptim=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);

rng default % For reproducibility
%options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata,'lb',lb,'ub',ub);%'options',options);
ms=MultiStart;
poptim=run(ms,problem,20);

f=poptim;
end

function f=sim2fit(NNsectorAges,params,v7,data64)
%tvec=[params(2),82,83,84];%80=last day of interest
%tvec=[params,78,79,80];%82,83,84];%80=last day of interest
tvec=[-44.3012,78,79,80];%82,83,84];%80=last day of interest
[pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNsectorAges,data64);%,params(1));
[simu,~]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,ones(63,1),tvec,0,data64,params);
%tout=simu(v7,1);
f=simu;%(v7);%Cumulative deaths
end