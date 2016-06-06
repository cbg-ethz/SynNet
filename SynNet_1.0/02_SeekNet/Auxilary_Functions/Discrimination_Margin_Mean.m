% This function reports the log10 fold change between the median of the "positive"
% and the median of the "negative" samples in the AnnoationCnt.
% Pej 2014

function Margin = Discrimination_Margin_Mean(AnnoationCnt, TrueAnnotation)
lAnnoationCnt = log10(AnnoationCnt);
m1 = mean(lAnnoationCnt( TrueAnnotation));
M0 = mean(lAnnoationCnt(~TrueAnnotation));

Margin = (m1) - (M0);
end