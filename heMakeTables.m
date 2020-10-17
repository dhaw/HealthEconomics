function f=heMakeTables(x1,x2,x3,x4,x1s,x2s,x3s,x4s,ddata)
nMonths=6;
nSchools=55;
x1=reshape(x1,63,nMonths);
x2=reshape(x2,63,nMonths);
x3=reshape(x3,63,nMonths);
x4=reshape(x4,63,nMonths);
x1s=reshape(x1s,63,nMonths);
x2s=reshape(x2s,63,nMonths);
x3s=reshape(x3s,63,nMonths);
x4s=reshape(x4s,63,nMonths);
xmin=ddata.xmin';
xmins=xmin;
xmin(nSchools)=.8;

T1=array2table([xmin,x1,x2,x3,x4]);
T1.Properties.VariableNames={'xmin','A1m1','A1m2','A1m3','A1m4','A1m5','A1m6','A2m1','A2m2','A2m3','A2m4','A2m5','A2m6','A3m1','A3m2','A3m3','A3m4','A3m5','A3m6','A4m1','A4m2','A4m3','A4m4','A4m5','A4m6'};
writetable(T1,'scenarioA.csv')

T2=array2table([xmins,x1s,x2s,x3s,x4s]);
T2.Properties.VariableNames={'xmin','B1m1','B1m2','B1m3','B1m4','B1m5','B1m6','B2m1','B2m2','B2m3','B2m4','B2m5','B2m6','B3m1','B3m2','B3m3','B3m4','B3m5','B3m6','B4m1','B4m2','B4m3','B4m4','B4m5','B4m6'};
writetable(T2,'scenarioB.csv')