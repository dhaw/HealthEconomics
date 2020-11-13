function f=makeTableXoptim(x)
scenLetter='A';
Hmax=12000;

lx=length(x);
numSect=63;
numInt=lx/numSect;

names=cellstr(num2str((1:numInt)'));

x=reshape(x,63,numInt);
Y=[(1:numSect)',x];
Y=array2table(Y);
Y.Properties.VariableNames=['Sector';names]';

writetable(Y,strcat('Scen',scenLetter,'Hmax',num2str(Hmax),'.csv'));