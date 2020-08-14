%X=[f1p1(:,1),sum(f1p1(:,2:end),2),sum(f1p4(:,2:end),2),sum(f1p7(:,2:end),2)];

X=[f1p1(:,1),f1p1(:,2:end),f1p4(:,2:end),f1p7(:,2:end)];
X(:,5)=X(:,2)+X(:,5);
X(:,10)=X(:,7)+X(:,10);
X(:,15)=X(:,12)+X(:,15);
X(:,[2,7,12])=[];

X(:,1)=discretize(X(:,1),62:7:X(end,1));
X(isnan(X(:,1)),:)=[];

%Y=[accumarray(X(:,1),X(:,2)),accumarray(X(:,1),X(:,3)),accumarray(X(:,1),X(:,4))];

Y=[accumarray(X(:,1),X(:,2)),accumarray(X(:,1),X(:,3)),accumarray(X(:,1),X(:,4)),accumarray(X(:,1),X(:,5)),accumarray(X(:,1),X(:,6)),accumarray(X(:,1),X(:,7)),accumarray(X(:,1),X(:,8)),accumarray(X(:,1),X(:,9)),accumarray(X(:,1),X(:,10)),accumarray(X(:,1),X(:,11)),accumarray(X(:,1),X(:,12)),accumarray(X(:,1),X(:,13))];

Y1=[p(:,1)+p(:,2)+.1114*p(:,3),(1-.1114)*p(:,3),p(:,4)];
Y2=[p(:,5)+p(:,6)+.1114*p(:,7),(1-.1114)*p(:,7),p(:,8)];
Y3=[p(:,9)+p(:,10)+.1114*p(:,11),(1-.1114)*p(:,11),p(:,12)];

figure; plot(Y)
Y=[(1:size(Y,1))',Y];

Y=[(1:size(Y1,1))',Y1,Y2,Y3];

Y=array2table(Y);

%Y.Properties.VariableNames={'Week','R<1.1','R<1.4','R<1.4','R<1.7'};

Y.Properties.VariableNames={'Week','R<1.1, a1','R<1.1, a2','R<1.1, a3','R<1.1, a4','R<1.4, a1','R<1.4, a2','R<1.4, a3','R<1.4, a4','R<1.7, a1','R<1.7, a2','R<1.7, a3','R<1.7, a4'};
Y.Properties.VariableNames={'Week','R<1.1, a1','R<1.1, a2','R<1.1, a3','R<1.4, a1','R<1.4, a2','R<1.4, a3','R<1.7, a1','R<1.7, a2','R<1.7, a3'};

writetable(Y,'projections.csv');