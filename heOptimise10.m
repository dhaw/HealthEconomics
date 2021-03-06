function [xoptim,yoptim,exitflag]=heOptimise10%(xoptim)%(thresh2)
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
%beta=0.0877;%beta64
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
%{
%Econ constrained by month:
kkron=eye(numInt-1);
kkron(end+1,end+1)=1/lastInt;%Account for last interval being longer
Gall=sparse(kron(kkron,G));
lx=size(Gall,1);
%Gvec=diag(Gall);
Z2=[Gall;-Gall];%<0,b
b2=[ball/numInt;zeros(lx,1)];
%
lb=zeros(lx,1);
ub=ones(lx,1);
X0=.1*rand(lx,1);
%}
%
%Econ constrained over whole period:
Z2=[repmat(G,1,numInt);repmat(-G,1,numInt)];
b2=[b';zeros(numSect,1)]+800;
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
%{
%options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e6,'algorithm','interior-point');%'sqp' %'interior-point'
[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub);%,nonlcon,options);
%}
%{
%rng default % For reproducibility
gs=GlobalSearch;
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon);
[xoptim,yoptim,exitflag]=run(gs,problem);
%}
%{
rng default % For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6,'algorithm','interior-point','UseParallel',true);%best death timeseries from public ONS
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
ms=MultiStart;
[xoptim,yoptim]=run(ms,problem,12);
exitflag='Not relevant for ms';
%}
%{ult
options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'UseVectorized', false);
[xoptim,yoptim,exitflag]=patternsearch(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
%}
%{
rng default
options=optimoptions('ga','UseParallel',true,'InitialPopulationRange',[lb'; ub']);
%opts.InitialPopulationRange=[lb'; X0'];%2 rows, nvar columns - lb;ub
[xoptim,yoptim,exitflag]=ga(fun1,378,Z2,b2,[],[],lb,ub,nonlcon,options);
%}
save('patt1xoutR1.mat','xoptim','yoptim','exitflag')
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

%{
G=[2047	-1455	0	-52	0	0	-5	0	-7	-5;
-701	27786	-4789	-8256	-1337	-715	-455	-2178	-6262	-1218;
-29	-651	16336	-426	-15	-399	-1579	-294	-602	-29;
-116	-1978	-112	47726	-475	-1468	-57	-1542	-1271	-188;
-39	-689	-220	-996	13742	-1732	-158	-952	-927	-349;
-116	-1513	-376	-1089	-239	18088	-2333	-758	-860	-124;
0	-156	-33	-1487	-299	-338	28938	-488	-716	-124;
-153	-3495	-1385	-4594	-2060	-3835	-699	29786	-3199	-1197;
-2	-444	-99	-540	-57	-212	-260	-490	41094	-23;
-6	-26	0	-16	-199	-135	-2	-64	-440	7754];
b=[1570,5628,36931,121561,23040,32041,75895,27509,116901,20598];
%}