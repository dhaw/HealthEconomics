function f=fitt0(NNsectorAges,ydata,datax)
%Data=matrix - vector of weekly deaths by week of year
xdata=(80:213)';%(1:wend)';
x0=-49.7725;%pmod,t0 %0
lb=-200;%50
ub=82;%92
fun=@(params,xdata)sim2fit(NNsectorAges,params,datax,xdata);%,v7
rng default % For reproducibility
%options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata,'lb',lb,'ub',ub);%'options',options);
ms=MultiStart;
poptim=run(ms,problem,100);
f=poptim;
end

function f=sim2fit(NNsectorAges,params,datax,xdata)%,v7
tvec=[params,86.3881,214,215];
[pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNsectorAges,datax);
[simu,~]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,datax.xmin,tvec,0,datax);
if simu(1,2)>1
    simu=[(1:simu-1)';simu];
else
    simu=simu(:,2);
end
f=simu(xdata);%(v7);%Cumulative deaths
end