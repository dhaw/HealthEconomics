function f=heMakeOneTable(x1,ddata,eecon,delta)%(x1,x2,x3,x4,x1s,x2s,x3s,x4s,ddata,eecon)
nMonths=3;
period=6/nMonths;

scen='B2delta0p';
xmin=ddata.xmin';
xfull=max(1,xmin);%Upper bound
xfull=repmat(xfull,nMonths,1);
%%
x1=reshape(x1-xfull,63,nMonths);
obj=repmat(eecon.obj,1,nMonths)*period;
x1=x1.*obj;
%%
T1=array2table(x1);
%T1.Properties.VariableNames={'xmin','A1m1','A1m2','A1m3','A1m4','A1m5','A1m6','A2m1','A2m2','A2m3','A2m4','A2m5','A2m6','A3m1','A3m2','A3m3','A3m4','A3m5','A3m6','A4m1','A4m2','A4m3','A4m4','A4m5','A4m6'};
T1.Properties.VariableNames={'HmaxPeriod1','HmaxPeriod2','HmaxPeriod3'};
%writetable(T1,strcat('scenario',scen,num2str(delta),'.csv'))
writetable(T1,'scenB2schoolsLessSus.csv')