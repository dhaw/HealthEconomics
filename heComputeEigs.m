function f=heComputeEigs(pr,beta,D,NNbar,ntot,S)%ntot=nbar
Deff=D.*repmat(S,1,ntot)./repmat(NNbar',ntot,1);
%
onesn=ones(ntot,1);
F=zeros(4*ntot,4*ntot);
F(1:ntot,ntot+1:end)=[pr.red*Deff,repmat(Deff,1,2)];
%vvec=kron([pr.sigma;pr.g1;pr.g2;pr.gX],ones(ntot,1));
%vvec(end-ntot+1:end)=vvec(end-ntot+1:end)+pr.h;
vvec=[pr.sigma*onesn;pr.g1*onesn;pr.g2*onesn;pr.gX+pr.h];
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=diag(-pr.sigma*(1-pr.p1)*onesn);
%V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*(1-pr.p2)*onesn);
V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*(1-pr.p2));
%V(3*ntot+1:4*ntot,1:ntot)=diag(-pr.sigma*pr.p1*pr.p2*onesn);
V(3*ntot+1:4*ntot,1:ntot)=diag(-pr.sigma*pr.p1*pr.p2);
GD=F/V;
f=eigs(GD);
f=beta*f(1);