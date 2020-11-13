function g=heSingleSim(xoptim)
load('NNs64.mat','NNs64')
load('ddata64.mat','ddata64')
%load('eecon64.mat','eecon64')
[pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNs64,ddata64);
lx=length(xoptim);
if lx==189
    tvec=[-49.7725   86.3881  245 306  367  426];
elseif lx==178
    tvec=[-49.7725   86.3881  245 275  306  336  367  398  426];
else
    error('Input not consistent with 63 sectors and 3x2-months or 6x1 month');
end
[f,g]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,ddata64);
plotMultiOutd(f,xoptim,tvec,ddata64);