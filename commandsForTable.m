function f=commandsForTable(frt,Rtmax)
ld2date='1stFeb';
%64 sector, one projection:
X=frt;%,frt1p1(:,2:end),frt1p25(:,2:end),frt1p4(:,2:end)];
X(:,67)=sum(X(:,2:64),2)+X(:,67);
X(:,2:64)=[];

X(:,1)=discretize(X(:,1),62:7:X(end,1)+7);
X(isnan(X(:,1)),:)=[];


Y=[accumarray(X(:,1),X(:,2)),accumarray(X(:,1),X(:,3)),accumarray(X(:,1),X(:,4)),accumarray(X(:,1),X(:,5))];
p=Y;
Y=[p(:,1)+p(:,2)+.1114*p(:,3),(1-.1114)*p(:,3),p(:,4)];
%figure; plot(Y)
%Y=[(1:size(Y,1))',Y];
factor=55.98/66.27;
Y=[(1:size(Y,1))',factor*Y];
f=sum(Y(:,2:end),2);

Y=array2table(Y);
%setHeading(Y,strcat('Hospital incidence, Rt(max)=',num2str(Rtmax)));

%Y.Properties.VariableNames={'Week','R<1.1','R<1.4','R<1.4','R<1.7'};
%Y.Properties.VariableNames={'Week','R<1.1, a1','R<1.1, a2','R<1.1, a3','R<1.1, a4','R<1.4, a1','R<1.4, a2','R<1.4, a3','R<1.4, a4','R<1.7, a1','R<1.7, a2','R<1.7, a3','R<1.7, a4'};
%Y.Properties.VariableNames={'Week','R<1.1, a1','R<1.1, a2','R<1.1, a3','R<1.25, a1','R<1.25, a2','R<1.25, a3','R<1.4, a1','R<1.4, a2','R<1.4, a3'};
Y.Properties.VariableNames={'Week','0-24','25-64','65+'};

writetable(Y,strcat('hosp_Rt',num2str(Rtmax),'LD',ld2date,'.csv'));
%{
%[f,g]=heRunCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta,repmat(ddata64.xmin',2,1),[tvec([1,2,3,6]),670],1,ddata64,.6345)
%Command to plot:
figure;
plot(X,'-','linewidth',2)
xlabel('Week')
ylabel('New hosp.')
axis tight
grid on
grid minor
box on
%legend('1','1.05')
end
%}
%%
%{
%X=[f1p1(:,1),sum(f1p1(:,2:end),2),sum(f1p4(:,2:end),2),sum(f1p7(:,2:end),2)];

X=[frt1p1(:,1),frt1p1(:,2:end),frt1p25(:,2:end),frt1p4(:,2:end)];
X(:,5)=X(:,2)+X(:,5);
X(:,10)=X(:,7)+X(:,10);
X(:,15)=X(:,12)+X(:,15);
X(:,[2,7,12])=[];

%10 sector:
X=[frt1p1(:,1),frt1p1(:,2:end),frt1p25(:,2:end),frt1p4(:,2:end)];
X(:,14)=sum(X(:,2:11),2)+X(:,14);
X(:,28)=sum(X(:,16:25),2)+X(:,28);
X(:,42)=sum(X(:,30:39),2)+X(:,42);
X(:,[2:11,16:25,30:39])=[];

X(:,1)=discretize(X(:,1),62:7:X(end,1));
X(isnan(X(:,1)),:)=[];

%Y=[accumarray(X(:,1),X(:,2)),accumarray(X(:,1),X(:,3)),accumarray(X(:,1),X(:,4))];

Y=[accumarray(X(:,1),X(:,2)),accumarray(X(:,1),X(:,3)),accumarray(X(:,1),X(:,4)),accumarray(X(:,1),X(:,5)),accumarray(X(:,1),X(:,6)),accumarray(X(:,1),X(:,7)),accumarray(X(:,1),X(:,8)),accumarray(X(:,1),X(:,9)),accumarray(X(:,1),X(:,10)),accumarray(X(:,1),X(:,11)),accumarray(X(:,1),X(:,12)),accumarray(X(:,1),X(:,13))];
p=Y;

Y1=[p(:,1)+p(:,2)+.1114*p(:,3),(1-.1114)*p(:,3),p(:,4)];
Y2=[p(:,5)+p(:,6)+.1114*p(:,7),(1-.1114)*p(:,7),p(:,8)];
Y3=[p(:,9)+p(:,10)+.1114*p(:,11),(1-.1114)*p(:,11),p(:,12)];

figure; plot(Y)
Y=[(1:size(Y,1))',Y];

factor=55.98/66.27;
Y=[(1:size(Y1,1))',factor*[Y1,Y2,Y3]];

Y=array2table(Y);

%Y.Properties.VariableNames={'Week','R<1.1','R<1.4','R<1.4','R<1.7'};
%Y.Properties.VariableNames={'Week','R<1.1, a1','R<1.1, a2','R<1.1, a3','R<1.1, a4','R<1.4, a1','R<1.4, a2','R<1.4, a3','R<1.4, a4','R<1.7, a1','R<1.7, a2','R<1.7, a3','R<1.7, a4'};
%Y.Properties.VariableNames={'Week','R<1.1, a1','R<1.1, a2','R<1.1, a3','R<1.25, a1','R<1.25, a2','R<1.25, a3','R<1.4, a1','R<1.4, a2','R<1.4, a3'};
Y.Properties.VariableNames={'Week','R<R_{LD}, a1','R<_{LD}, a2','R<_{LD}, a3','R<1.05, a1','R<1.05, a2','R<1.05, a3','R<1.1, a1','R<1.1, a2','R<1.1, a3'};

writetable(Y,'projections.csv');
%}