function f=heMakeDs(NN,x,datax,int)%NNfrac,x4)
%NN - vector of populations including non-working. 
%x - proportion of each sector open - not including non-working.
lx=length(x);%Number of sectors
%%
%Size of C explicit in here
if int==0
    background=datax.comm(1);
elseif int==1
    background=datax.comm(1);
else
    background=datax.comm(1);
end
C=[1.7112    0.4927    0.4049    0.0293;
    0.3919    2.5527    0.4753    0.0348;
    0.3234    0.5548    0.8996    0.0728;
    0.0528    0.1904    0.3744    0.3830];
%C=C*data64.comm(1)/mean(mean(C));%Match community transmission

adInd=3;%Adult index
%
lc=length(C);
C=C*background/sum(C(3,:));
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
elseif length(NN)==67
    ln=length(NN);
    NNrel=NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));
    %Make A:
    matA=zeros(ln,ln);
    matA(lx+1:end,lx+1:end)=C;
    matA(1:lx,lx+1:end)=repmat(CworkRow,lx,1);
    matA(:,[1:lx,lx+adInd])=repmat(matA(:,lx+adInd),1,lx+1).*repmat(NNrel',ln,1);
    %Make B and C:
    valB=datax.B;
    valC=datax.C;
    if int>0
        valB=valB.*(1-datax.wfhAv);
        valC=valC.*(1-datax.wfhAv);
    end
    valB(lx+1:ln)=0;
    valC(lx:1:ln)=0;
    %
    x(lx+1:lx+lc)=0;
    matB=diag(x.*valB');
    matB(ln,ln)=0;
    NNrep=repmat(NN'/sum(NN),lx+lc,1);
    matC=repmat(x.*valC',1,ln).*NNrep;
    
    %Modify depending on x:
    %Education:
    matA(lx+1,lx+1)=matA(lx+1,lx+1)+datax.schoolA1*x(55);
    matA(lx+2,lx+2)=matA(lx+2,lx+2)+datax.schoolA2*x(55);
    %Hospitality:
    sects=[58,59,60,62];
    psub=datax.NNsector(sects)';
    psub=sum(psub.*x(sects))/sum(psub);
    matA([1:lx,lx+3:ln],:)=matA([1:lx,lx+3:ln],:)+NNrep([1:lx,lx+3:ln],:)*psub*datax.hospA34(1);
    matA(lx+2,:)=matA(lx+2,:)+NNrep(lx+2,:)*psub*datax.hospA2;%(1)
    %Transport:
    if int==0
        matA(1:lx,1:lx)=matA(1:lx,1:lx)+NNrep(1:lx,1:lx)*datax.travelA3(1);
    else
        matA(1:lx,1:lx)=matA(1:lx,1:lx)+NNrep(1:lx,1:lx)*datax.travelA3(1).*repmat(1-datax.wfhAv,lx,1).*repmat(1-datax.wfhAv',1,lx);
    end
    %
    f=matA+matB+matC;
elseif length(NN)==5%Single sector
    ln=length(NN);
    NNrel=NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));
    %Make A:
    matA=zeros(ln,ln);
    matA(lx+1:end,lx+1:end)=C;
    matA(1:lx,lx+1:end)=repmat(CworkRow,lx,1);
    matA(:,[1:lx,lx+adInd])=repmat(matA(:,lx+adInd),1,lx+1).*repmat(NNrel',ln,1);
    %Make B and C:
    valB=datax.B;
    valC=datax.C;
    if int>0
        valB=valB.*(1-datax.wfhAv);
        valC=valC.*(1-datax.wfhAv);
    end
    valB(lx+1:ln)=0;
    valC(lx:1:ln)=0;
    %
    x(lx+1:lx+lc,1)=0;
    matB=diag(x.*valB');
    matB(ln,ln)=0;
    NNrep=repmat(NN'/sum(NN),lx+lc,1);
    matC=repmat(x.*valC',1,ln).*NNrep;
    
    %Modify depending on x:
    %Education:
    matA(lx+1,lx+1)=matA(lx+1,lx+1)+datax.propschools*datax.schoolA1*x(1);
    matA(lx+2,lx+2)=matA(lx+2,lx+2)+datax.propschools*datax.schoolA2*x(1);
    %Hospitality:
    matA([1:lx,lx+3:ln],:)=matA([1:lx,lx+3:ln],:)+NNrep([1:lx,lx+3:ln],:)*datax.prophosp*datax.hospA34;
    matA(lx+2,:)=matA(lx+2,:)+NNrep(lx+2,:)*datax.prophosp*datax.hospA2;%(1)
    %Transport:
    if int==0
        matA(1:lx,1:lx)=matA(1:lx,1:lx)+NNrep(1:lx,1:lx)*datax.travelA3(1);
    else
        matA(1:lx,1:lx)=matA(1:lx,1:lx)+NNrep(1:lx,1:lx)*datax.travelA3(1).*repmat(1-datax.wfhAv,lx,1).*repmat(1-datax.wfhAv',1,lx);
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