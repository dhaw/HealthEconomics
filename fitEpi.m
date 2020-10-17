function f=fitEpi(NNsectorAges,ydata,datax)

NNsectorAges=NNsectorAges*56286961/66796807;
datax.NNsector=datax.NNsector*56286961/66796807;

%Data=matrix - vector of weekly deaths by week of year
%hePrepCovid19 - feed in R0
%heSimCovid19 - output deaths
%wend=11;
xdata=(80:213);%20th March to 30th Jun;%+60;%(41:152)';%(1:wend)';
%ydata=cumsum(ydata(1:wend));
%v7=(5:7:wend*7-2)';
%
x0=[.3,0,83,3];%pmod,t0 %0
lb=[0,-100,0,2.5];%50
ub=[1,91,92,3.5];%92
%}
%{
x0=[.3,0,83];%pmod,t0 %0
lb=[0,-100,0];%50
ub=[1,91,92];%92
%}
%{
x0=.4;
lb=0;%-100;
ub=1;%100;
%}
fun=@(params,xdata)sim2fit(NNsectorAges,params,datax,xdata);%,v7
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

function f=sim2fit(NNsectorAges,params,datax,xdata)%,v7
tvec=[params(2),params(3),214,215];%78,79,80];%80=last day of interest params(3) -> 83  122,123
%tvec=[-49.7725,86.3881,153,154];
%tvec=[params,78,79,80];%82,83,84];%80=last day of interest
%tvec=[-44.3012,78,79,80];%82,83,84];%80=last day of interest
[pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNsectorAges,datax,params(4));
[simu,~]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,datax.xmin,tvec,0,datax,params(1));
%tout=simu(v7,1);
if simu(1,1)>1%(1,2)****?
    %simu=[(1:simu-1)';simu(:,2)];
    simu=[zeros(simu(1,1)-1,1);simu(:,2)];
else
    simu=simu(:,2);
end
f=simu(xdata);%(v7);%Cumulative deaths
end