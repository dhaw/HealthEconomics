function f=he64to1(NNs64,data64)
data1=struct;
NNwork=sum(data64.NNsector);
props=data64.NNsector/NNwork;
%Need A,B,ctt,WFH,furl,pre

NNsum=sum(data64.NNsector);
data1.B=sum(data64.B.*props);
data1.C=sum(data64.C.*props);
data1.Bdur=sum(data64.Bdur.*props);
data1.Cdur=sum(data64.Cdur.*props);
data1.contAv=sum(data64.contAv.*props);
data1.wfhAv=sum(data64.wfhAv.*props);
data1.furlAv=sum(data64.furlAv.*props);
data1.xmin=sum(data64.xmin.*props);

data1.comm=data64.comm(1);

NNs64sum=sum(NNs64);
propA2=NNs64(65)/NNs64sum;
propA3=sum(NNs64([1:63,66]))/NNs64sum;
propA4=NNs64(67)/NNs64sum;
data1.travelA3=data64.travelA3(1);%*propA3;
%data1.hosp=data64.hospA2*propA2+data64.hospA34*(propA3+propA4);
data1.hospA2=data64.hospA2;
data1.hospA34=data64.hospA34(1);

data1.schoolA1=data64.schoolA1;
data1.schoolA2=data64.schoolA2;
data1.prophosp=sum(data64.NNsector([58,59,60,62]))/NNsum;
data1.propschools=data64.NNsector(55)/NNsum;

f=data1;