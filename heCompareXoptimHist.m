function f=heCompareXoptimHist(x1,x2,x3,x4,ddata)%(x1,x1s,x2,x2s,x3,x3s,x4,x4s)%x2 schools forced open
NN=ddata.NNsector;
bins=[0,cumsum(NN')];
mids=bins(1:end-1)+.5+diff(bins);

%hmax=1;
H={'18,057','28,115','38,172','53,365'};

inds=repelem((1:10)',[3,23,1,9,4,3,1,9,4,6]');
inds=repmat(inds,1,4);

topn=5;
lx=length(x1)/6;
X=zeros(lx,4);
Y=zeros(topn,4);

x1=reshape(x1,lx,6);
X1=mean(x1,2).*NN;

x2=reshape(x2,lx,6);
X2=mean(x2,2).*NN;

x3=reshape(x3,lx,6);
X3=mean(x3,2).*NN;

x4=reshape(x4,lx,6);
X4=mean(x4,2).*NN;


%numInds=5;
%inds=xdiffSum(1:numInds,1);
%indSch=55;%find(xdiffSum(:,1)==55);


fs=10; lw=2; ms=8;
col1=.5*[1,1,1];
%cmap=.8*flipud(winter(numInds));
cmap=jet(10);
maxy=max(max(X))/6;
miny=min(min(X))/6;

figure

subfigure(2,2,1)
X=X1;
hold on
for i=1:10
    histogram(X.*mids,bins)
end
xlabel('H_{max}')
ylabel('\Delta x_i')
set(gca,'fontsize',fs)
axis([0,nvec(end)+1,0,maxy])
grid on
grid minor
box on

subfigure(2,2,2)
X=X2;
hold on
for i=1:10
    histogram(X.*mids,bins)
end
xlabel('H_{max}')
ylabel('\Delta x_i')
set(gca,'fontsize',fs)
axis([0,nvec(end)+1,0,maxy])
grid on
grid minor
box on

subfigure(2,2,3)
X=X3;
hold on
for i=1:10
    histogram(X.*mids,bins)
end
xlabel('H_{max}')
ylabel('\Delta x_i')
set(gca,'fontsize',fs)
axis([0,nvec(end)+1,0,maxy])
grid on
grid minor
box on

subfigure(2,2,4)
X=X4;
hold on
for i=1:10
    histogram(X.*mids,bins)
end
xlabel('H_{max}')
ylabel('\Delta x_i')
set(gca,'fontsize',fs)
axis([0,nvec(end)+1,0,maxy])
grid on
grid minor
box on