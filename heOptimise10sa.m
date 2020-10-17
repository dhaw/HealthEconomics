function [xoptim,yoptim,exitflag]=heOptimise10sa%(xoptim)%(thresh2)
%%
load('NNs10.mat','NNs10')
NNsector=NNs10;
load('ddata10.mat','ddata10')
%load('GAxoutR1schools.mat.mat','xoptim')
%%
hospThresh=[18057,1];%38172	53365	74966	82966 % 18057   28115   38172	53365

t0=-59.4613;%Re-fitted as single parameter %-49.7725;%-53.837;
t1=86.3881;%87.1771;
tend=720;
months=[1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,tend];
tvec=[t0,t1,months(9:15)];%Fit tvec(1)=t0 %11
lastInt=1;%(tend-tvec(end-1))/30;
numInt=length(tvec)-3;
numAges=4;
numSect=length(NNsector)-numAges;
%%
%Tune epi-midel to pre-lockdown:
[pr,NN,n,nbar,na,NNbar,NNrep,Din,beta]=hePrepCovid19(NNsector,ddata10);
%%
%
G=[1962	-1022	-5	-30	0	0	-8	-5	-4	-3;
-430	38558	-3199	-4783	-496	-358	-316	-980	-2944	-428;
-21	-473	16041	-343	0	-303	-1817	-160	-508	-15;
-256	-5406	-1125	43670	-494	-1192	-160	-1725	-2230	-251;
-31	-419	-213	-853	14186	-1390	-252	-1019	-720	-231;
-96	-1298	-338	-979	-247	18587	-2051	-665	-695	-106;
0	-74	-16	-1530	-263	-294	28354	-451	-649	-102;
-104	-1954	-1124	-3699	-1660	-2216	-762	27555	-2421	-710;
-2	-140	-108	-248	-72	-187	-347	-482	39791	-14;
-2	-16	-1	-6	-138	-107	-7	-57	-351	7210];
b=[5318,147744,74405,184988,54354,72676,149844,77424,229146,39148];
objFun=[922,18515,8925,27340,9296,9765,22150,20185,26164,4914];%Monthly
%}
ball=repmat(b',numInt,1);
objFun=repmat(objFun',numInt,1);
%%
%
%Econ constrained over whole period:
Z2=[repmat(G,1,numInt);repmat(-G,1,numInt)];
b2=[b';zeros(numSect,1)];
lx=numInt*numSect;
lb=repmat(ddata10.xmin',numInt,1);%zeros(lx,1);
ub=ones(lx,1);
%X0=repmat(ddata10.xmin',numInt,1);%zeros(lx,1);
X0=[repmat(ones(10,1),5,1);repmat(ddata10.xmin',1,1)];
%{
%Schools open fully in last 2 months:
xlb=repmat(ddata10.xmin',numInt,1);
%newx=(1-ddata10.propschools)*ddata10.xmin(9)+ddata10.propschools*.8;
newx=(1-ddata10.propschools)*0.7991+ddata10.propschools*.8;%Weighted average of non-education=0.7991
xlb((0:5)*10+9,1)=newx;
%xlb(5*10+9,1)=newx;
lb=xlb;
%X0=xlb;%repmat(xlb',numInt,1);%zeros(lx,1);
X0((0:5)*10+9,1)=xnew;
%}
%X0=xoptim;
%xmin6=repmat(ddata10.xmin',6,1);
%X0=xmin6+.99*(xoptim-xmin6);
%Schools open:
%X0(51:60)=xmin6(51:60);
%}
%%
fun1=@(Xit)econGDP(objFun,Xit);
%fun2=@(Xit)econConstraint(Z,Xit);
nonlcon=@(Xit)epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata10);
fun2=@(Xit)forSimAnn(objFun,Z2,b2,pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata10);
rng default % For reproducibility
%options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'algorithm','interior-point','UseParallel',true);
[xoptim,yoptim,exitflag,output]=simulannealbnd(fun2,X0,lb,ub);%,options);
%%
save('SA1xoutR1.mat','xoptim','yoptim','exitflag','output')
%}
%%

end

function f=econGDP(G,Xit)%,lx)
f=-sum(G.*Xit);
%f=-sum(sum(G*Xit));
%{
y=0;
for i=1:lx[
    y=y+sum(G*X(:,i));
end
f=-y;
%}
end

function [c,cex]=epiConstraint(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,data10)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0,data10);
c=max(h-hospThresh);
cex=[];
end

function f=forSimAnn(obj,A,b,pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,ddata10)
[~,h]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,0,ddata10);
if max(h-hospThresh)>0
    f=0;
elseif max(A*Xit-b)>0
    f=0;
else
    f=-sum(obj.*Xit);
end
end