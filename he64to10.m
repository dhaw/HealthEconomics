function [NNs10,data10]=he64to10(NNs64,data64,agg)
%Everything a row vector - except agg
data10=struct;
data10.NNsector=accumarray(agg,data64.NNsector')';
[freq,~]=groupcounts(agg);
tots=repelem(data10.NNsector',freq);
props=data64.NNsector'./tots;
NNs10=[data10.NNsector';NNs64(64:end)];

data10.xmin=accumarray(agg,data64.xmin'.*props)';

%Need A,B,ctt,WFH,furl,pre

data10.B=accumarray(agg,data64.B'.*props)';
data10.C=accumarray(agg,data64.C'.*props)';
%data10.Bdur=accumarray(agg,data64.Bdur'.*props)';
%data10.Cdur=accumarray(agg,data64.Cdur'.*props)';
data10.contAv=accumarray(agg,data64.contAv'.*props)';
data10.wfhAv=accumarray(agg,data64.wfhAv'.*props)';
data10.furlAv=accumarray(agg,data64.furlAv'.*props)';
data10.xmin=accumarray(agg,data64.xmin'.*props)';

data10.comm=data64.comm(1);

NNs10sum=sum(NNs10);
propA2=NNs64(12)/NNs10sum;%Proportion school age
propA3=sum(NNs10([1:10,13]))/NNs10sum;%Proportion working age
propA4=NNs64(14)/NNs10sum;%Proportion old age
data10.travelA3=data64.travelA3(1);%*propA3;
%data1.hosp=data64.hospA2*propA2+data64.hospA34*(propA3+propA4);
data10.hospA2=data64.hospA2;
%data10.hospA34=data64.hospA34(1);
data10.hospA3=data64.hospA3(1);
data10.hospA4=data64.hospA4(1);

data10.schoolA1=data64.schoolA1;
data10.schoolA2=data64.schoolA2;
data10.prophosp=sum(data64.NNsector([58,59,60,62]))/data10.NNsector(10);
data10.propschools=data64.NNsector(55)/data10.NNsector(9);