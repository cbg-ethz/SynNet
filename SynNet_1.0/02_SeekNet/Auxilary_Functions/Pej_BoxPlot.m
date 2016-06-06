function Pej_BoxPlot(X)
if size(X,1)==1; X = X';end

L = size(X,1);
boxplot(X,'colors', 'k', 'symbol','')
hold on
for i = 1:size(X,2)
    tmpX= min(max(randn(L,1)/15,-.35),.35) + i;
    plot(tmpX,X(:,i), 's', 'color', [1 1 1]*.5, 'markerfacecolor', [1 1 1]*.5, 'markersize', 3)
end
end