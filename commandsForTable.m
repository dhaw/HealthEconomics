X=[f1p1(:,1),sum(f1p1(:,2:end),2),sum(f1p4(:,2:end),2),sum(f1p7(:,2:end),2)];
X(:,1)=discretize(X(:,1),62:7:X(end,1));
X(isnan(X(:,1)),:)=[];
Y=[accumarray(X(:,1),X(:,2)),accumarray(X(:,1),X(:,3)),accumarray(X(:,1),X(:,4))];
figure; plot(Y)
Y=[(1:size(Y,1))',Y];
Y=array2table(Y);
Y.Properties.VariableNames={'Week','R<1.1','R<1.4','R<1.7'};
writetable(Y,'projections.csv');