function f=fitEpi(NNsectorAges,ydata,datax)
%Data=matrix - vector of weekly deaths by week of year
%hePrepCovid19 - feed in R0
%heSimCovid19 - output deaths
%wend=11;
xdata=(41:121)';%(1:wend)';
%ydata=cumsum(ydata(1:wend));
%v7=(5:7:wend*7-2)';
%
x0=[.5,0,92,3.5];%pmod,t0
lb=[0,-100,0,2.5];
ub=[1,100,150,5];
%}
%{
x0=0;
lb=-100;%-100;
ub=100;%100;
%}
fun=@(params,xdata)sim2fit(NNsectorAges,params,datax);%,v7
%fun=@(params)sim2fit(NNsectorAges,params,v7,datax);
%poptim=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
%
rng default % For reproducibility
%options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata,'lb',lb,'ub',ub);%'options',options);
ms=MultiStart;
poptim=run(ms,problem,100);
%}
%{
fun=@(params)sim2fit(NNsectorAges,params,v7,datax);
rng default
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('fmincon','x0',x0,'objective',fun,'lb',lb,'ub',ub,'options',options);
ms=MultiStart;
[poptim,~]=run(ms,problem,20);
%}
f=poptim;
end

function f=sim2fit(NNsectorAges,params,datax)%,v7
tvec=[params(2),params(3),122,123];%78,79,80];%80=last day of interest
%tvec=[params,78,79,80];%82,83,84];%80=last day of interest
%tvec=[-44.3012,78,79,80];%82,83,84];%80=last day of interest
[pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNsectorAges,datax,params(4));
[simu,~]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,datax.xmin,tvec,0,datax,params(1));
%tout=simu(v7,1);
if simu(1,2)>1
    simu=[(1:simu-1)';simu];
else
    simu=simu(:,2);
end
f=simu(41:121);%(v7);%Cumulative deaths
end