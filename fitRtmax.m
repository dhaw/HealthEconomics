function f=fitRtmax(NNsectorAges,datax,RmaxIn)
%Data=matrix - vector of weekly deaths by week of year
%RmaxIn=1.4;
x0=.56;%pmod,t0 %0
lb=0;
ub=2;
xdata=1;
fun=@(params)sim2fit(NNsectorAges,params,datax);
fun2=@(params,xdata)nonlcon(NNsectorAges,params,datax,RmaxIn,xdata);
%rng default % For reproducibility
%options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);
%problem=createOptimProblem('fmincon','x0',x0,'objective',fun,'lb',lb,'ub',ub,'nonlcon',fun2);
%poptim=lsqcurvefit(fun2,x0,xdata,0,lb,ub);%,options);
problem=createOptimProblem('lsqcurvefit','objective',fun2,'x0',x0,'xdata',xdata,'ydata',0,'lb',lb,'ub',ub);
ms=MultiStart;
poptim=run(ms,problem,10);
f=poptim;
end

function f=sim2fit(NNsectorAges,params,datax)
tvec=[0.5262,-49.7725,245,670];%[-49.7725,86.3881,245,670];%[-59.4613,86.3881,122,670];
%1st March=day 61
%62+78*7=608
[pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNsectorAges,datax);
[Rmax,~]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,datax.xmin',tvec,0,datax,params);%pmod
f=-Rmax;
end

function f=nonlcon(NNsectorAges,params,datax,RmaxIn,xdata)%(params,RmaxIn)%
%x=fun(params);%sim2fit(NNsectorAges,params,datax);
tvec=[-49.7725,86.3881,245,336,670];%[-59.4613,86.3881,122,670];
[pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNsectorAges,datax);
[~,Rmax]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,repmat(datax.xmin',2,1),tvec,0,datax,params);%pmod
f=Rmax(xdata)-RmaxIn;
end