function WeightedPerf = Get_Weighted_Performance(Function_Output, Annotation, SortedEliteIdx)
global Consts
NF = size(Function_Output,1);

if nargin<3
    SortedEliteIdx = 1:NF;
    PushtoTop = false;
else
    PushtoTop = true;
end
NE = length(SortedEliteIdx);
nAnnotation = ~Annotation;

switch Consts.AnalysisMode
    case 'B'
        Clasification_pwr =  double(Function_Output == repmat(Annotation, NF, 1));
    case 'C'
        E_Pos = mean(Function_Output(:, Annotation),2);
        E_Neg = mean(Function_Output(:,nAnnotation),2);
        Clasification_pwr                = zeros(size(Function_Output));
        Clasification_pwr(:, Annotation) =  Function_Output(:, Annotation) - repmat(E_Neg, 1, sum( Annotation));
        Clasification_pwr(:,nAnnotation) = -Function_Output(:,nAnnotation) + repmat(E_Pos, 1, sum(nAnnotation));
end


% Re-weighting Samples
Elite_margin              = mean(Clasification_pwr(SortedEliteIdx,:),1);
Elite_margin              = Elite_margin - mean(Elite_margin);

SampleWeights = 10.^(-Elite_margin);
SampleWeights(nAnnotation) = SampleWeights(nAnnotation) * (sum(Annotation)/sum(nAnnotation));
SampleWeights = SampleWeights / mean(SampleWeights);
WeightedPerf  = max(Clasification_pwr,0) * SampleWeights'; % max(Clasification_pwr,0), because we just care about the postivie potential of the functions

if PushtoTop
    Offset = linspace(1,.1,length(SortedEliteIdx));
    WeightedPerf(SortedEliteIdx) = max(WeightedPerf)+Offset; % To keep the order as it was defined by the SortedEliteIdx
end
end