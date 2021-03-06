function f=heMakeDs64(NN,x)%NNfrac,x4)
%NN - vector of populations including non-working. 
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