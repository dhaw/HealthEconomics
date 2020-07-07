function f=heMakeDs(NN,x,data64,int)%NNfrac,x4)
%NN - vector of populations including non-working. 
%x - proportion of each sector open - not including non-working.
lx=length(x);%Number of sectors
%%
%Size of C explicit in here
background=data64.comm(1);
C=[1.7112    0.4927    0.4049    0.0293;
    0.3919    2.5527    0.4753    0.0348;
    0.3234    0.5548    0.8996    0.0728;
    0.0528    0.1904    0.3744    0.3830];
%C=C*data64.comm(1)/mean(mean(C));%Match community transmission

adInd=3;%Adult index
%
lc=length(C);
C=C*background/C(3,3);
CworkRow=C(3,:);
%%
%if length(NN)==lc+1%*1Dage
    %workContacts=27.4273;
    %f=[Cwork+x*workContacts,0;0,Cwork];
%else
if length(NN)==14
    ln=length(NN);
    NNrel=NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));
    %valA=Cwork;
    matA=zeros(lx+lc,lx+lc);
    matA(lx+1:end,lx+1:end)=C;
    matA(1:lx,lx+1:end)=repmat(CworkRow,lx,1);
    matA(:,[1:lx,lx+adInd])=repmat(matA(:,lx+adInd),1,lx+1).*repmat(NNrel',lx+lc,1);
    valB=[1.6,3.9,2.1,1.9,2.3,3,3,1.5,.9,2.3];
    valB(lx:1:lx+lc)=0;%"lx:1" -> "lx+1"
    valC=[4.9,17,58.2,86.3,21.6,13.8,13.8,17.7,46.8,21.6];
    valC(lx:1:lx+lc)=0;
    %NNrel=repmat(NNrel',lx+1,1);
    %matA=valA*NNrel;
    x(lx+1:lx+lc)=0;
    matB=diag(x.*valB');
    matB(ln,ln)=0;
    matC=repmat(x.*valC',1,ln).*repmat(NN'/sum(NN),lx+lc,1);
    %matC(ln,ln)=0;
    f=matA+matB+matC;
elseif length(NN)==64
    ln=length(NN);
    NNrel=NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));
    %Make A:
    matA=zeros(ln,ln);
    matA(lx+1:end,lx+1:end)=C;
    matA(1:lx,lx+1:end)=repmat(CworkRow,lx,1);
    matA(:,[1:lx,lx+adInd])=repmat(matA(:,lx+adInd),1,lx+1).*repmat(NNrel',ln,1);
    %Make B and C:
    valB=data64.B;
    valC=data64.C;
    valC(lx:1:lx+lc)=0;
    if int==1
        valB=valB.*(1-data64.wfhAv);
        valC=valC.*(1-data64.wfhAv);
    end
    valB(lx+1:ln)=0;
    valC(lx:1:lx+lc)=0;
    %
    x(lx+1:lx+lc)=0;
    matB=diag(x.*valB');
    matB(ln,ln)=0;
    matC=repmat(x.*valC',1,ln).*repmat(NN'/sum(NN),lx+lc,1);
    
    %Modify depending on x:
    %Education:
    matA(lx+1,lx+1)=matA(lx+1,lx+1)+data64.schoolA1*x(55);
    matA(lx+2,lx+2)=matA(lx+2,lx+2)+data64.schoolA2*x(55);
    %Hospitality:
    
    %Transport:
    if int==0
        matA(1:lx,1:lx)=matA(1:lx,1:lx)*data.travelA3;
    else
        matA(1:lx,1:lx)=matA(1:lx,1:lx)*data.travelA3*repmat(1-data64.wfhAv,lx,1).*repmat(1-data64.wfhAv',1,lx);
    end
    %
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