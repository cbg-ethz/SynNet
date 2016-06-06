function M = Pej_Mean_withNaNs(X)
M = mean(X(~isnan(X)));
end