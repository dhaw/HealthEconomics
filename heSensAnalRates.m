function H=heSensAnalRates(NNsx,ddatax,xoptim,tvec)
%
%heSim output: f=t,S,H
%tvec=[-53.8370   87.1771  122.0000  153.0000  183.0000  214.0000  245.0000  275.0000  306.0000];
sdfact=.05;
numIter=1000;
%
F=zeros(306,numIter);
%G=F;
H=F;
parfor i=1:numIter
    ddataxi=ddatax;
    %{
    ddataxi.B=normrnd(ddataxi.B,sdfact*ddataxi.B);
    ddataxi.C=normrnd(ddataxi.C,sdfact*ddataxi.C);
    %}{
    ddataxi.schoolA1=normrnd(ddataxi.schoolA1,sdfact*ddataxi.schoolA1);
    ddataxi.schoolA2=normrnd(ddataxi.schoolA2,sdfact*ddataxi.schoolA2);
    %}
    %{
    ddataxi.travelA3=normrnd(ddataxi.travelA3,sdfact*ddataxi.travelA3);
    ddataxi.hospA2=normrnd(ddataxi.hospA2,sdfact*ddataxi.hospA2);
    ddataxi.hospA3=normrnd(ddataxi.hospA3,sdfact*ddataxi.hospA3);
    ddataxi.hospA4=normrnd(ddataxi.hospA4,sdfact*ddataxi.hospA4);
    %}
    %
    ddataxi.comm=normrnd(ddataxi.comm,sdfact*ddataxi.comm);
    %}
    [pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(NNsx,ddataxi);
    [f1,~]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,ddataxi);
    %F(:,i)=f1(:,1);
    %G(:,i)=f1(:,2);
    H(:,i)=f1(:,3);
end