function f=heMakeTables2(x1,x2,x3,x4,x5,x6,x7,ddata,eecon)%(x1,x2,x3,x4,x1s,x2s,x3s,x4s,ddata,eecon)
nMonths=3;
period=6/nMonths;

obj=eecon.obj;
scen='All';
x0a=ddata.xmin';
x0b=x0a;
x0b(55)=.8;
xfull=max(1,ddata.xmin');

%%
lda=(x0a-xfull).*obj*6;
ldb=(x0b-xfull).*obj*6;

xfull=repmat(xfull,nMonths,1);
x1=reshape(x1-xfull,63,nMonths);
x1=sum(x1,2).*obj*period;
x2=reshape(x2-xfull,63,nMonths);
x2=sum(x2,2).*obj*period;
x3=reshape(x3-xfull,63,nMonths);
x3=sum(x3,2).*obj*period;
x4=reshape(x4-xfull,63,nMonths);
x4=sum(x4,2).*obj*period;
x5=reshape(x5-xfull,63,nMonths);
x5=sum(x5,2).*obj*period;
x6=reshape(x6-xfull,63,nMonths);
x6=sum(x6,2).*obj*period;

x7=reshape(x7-repmat(xfull,period,1),63,6);
x7=sum(x7,2).*obj;


%%
T=array2table([(1:63)',lda,ldb,x1,x2,x3,x4,x5,x6,x7]);
T.Properties.VariableNames={'Sector','LDA','LDB','A1','A2','A3','B1','B2','B3','B2(6x1)'};
writetable(T,'AllScenariosGVA.csv')