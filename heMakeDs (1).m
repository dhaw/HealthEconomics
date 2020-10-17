function f=heMakeDs(NN,x)%NNfrac,x4)
%NN - vector of populations including non-working. 
%x - proportion of each sector open - not including non-working. 
if length(NN)==2
    %{
    %workContacts=27.4273;
    %f=[7.2+x*workContacts,0;0,7.2];
    %
    valA=7.2;
    valB=.913;%[1.6,3.9,2.1,1.9,2.3,3,3,1.5,.9,2.3,0]';
    valC=22.0594;%[4.9,17,58.2,86.3,21.6,13.8,13.8,17.7,46.8,21.6,0]';
    %valB=sum(NN.*(valB))/sum(NN);
    %valC=sum(NN.*(valC))/sum(NN);
    propWork=NN(1)/sum(NN)*x;
    f=[propWork*(valA+valC)+x*valB,(1-propWork)*(valA+valC);0,valA];%propWork*valA,(1-propWork)*valA];
    %}
    %Pick a sector:
    sector=10;
    pw=NN(1)/sum(NN)*x;
    valA=7.2;
    valB=[1.6,3.9,2.1,1.9,2.3,3,3,1.5,.9,2.3,0]; valB=valB(sector);
    valC=[4.9,17,58.2,86.3,21.6,13.8,13.8,17.7,46.8,21.6,0]; valC=valC(sector);
    f=[pw*(valA+valC)+x*valB,(1-pw)*(valA+valC);0,valA];
else
    ln=length(NN);
    NNrel=NN/sum(NN);
    x(end+1)=0;

    valA=7.2;
    valB=[1.6,3.9,2.1,1.9,2.3,3,3,1.5,.9,2.3,0];
    valC=[4.9,17,58.2,86.3,21.6,13.8,13.8,17.7,46.8,21.6,0];

    NNrel=repmat(NNrel',ln,1);
    matA=valA*NNrel;
    matB=diag(x.*valB');
    matB(ln,ln)=0;
    matC=repmat(x.*valC',1,ln);
    matC(ln,ln)=0;

    f=matA+matB+matC;
end
%{
%Toy example:
baseRate=min(.1+x4,1);%x(4)=serveces
ln=length(NN);
NNsum=sum(NN);
NNrel=NN./NNsum;%Column
D=repmat(NNrel',ln,1)*baseRate;
D(1:ln-1,1:ln-1)=D(1:ln-1,1:ln-1)+diag(NNfrac);
f=D;
%}