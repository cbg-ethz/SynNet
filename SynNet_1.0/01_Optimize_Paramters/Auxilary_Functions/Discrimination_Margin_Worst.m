% This function reports the log2 fold change between the lowers "positive"
% and the highes "negative" sample in the AnnoationCnt.
% Pej 2014

function Margin = Discrimination_Margin_Worst(AnnoationCnt, TrueAnnotation)
m1 = min(AnnoationCnt( TrueAnnotation));
M0 = max(AnnoationCnt(~TrueAnnotation));

Margin = log2(m1) - log2(M0);
end