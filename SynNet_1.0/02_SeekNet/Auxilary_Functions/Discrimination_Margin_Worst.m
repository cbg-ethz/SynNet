% This function reports the log2 fold change between the lowers "positive"
% and the highes "negative" sample in the AnnoationCnt.
% Pej 2014

function Margin = Discrimination_Margin_Worst(AnnoationCnt, TrueAnnotation)
m1 = min(AnnoationCnt(:, TrueAnnotation),[],2);
M0 = max(AnnoationCnt(:,~TrueAnnotation),[],2);
Margin = log10(m1) - log10(M0);
end