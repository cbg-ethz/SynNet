% This function reports the log2 fold change between the median of the "positive"
% and the median of the "negative" samples in the AnnoationCnt.
% Pej 2014

function Margin = Discrimination_Margin_Median(AnnoationCnt, TrueAnnotation)
m1 = median(AnnoationCnt( TrueAnnotation));
M0 = median(AnnoationCnt(~TrueAnnotation));

Margin = log10(m1) - log10(M0);
end