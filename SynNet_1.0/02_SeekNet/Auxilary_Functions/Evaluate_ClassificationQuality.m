function [AverageMargin, WeightedPerf] = Evaluate_ClassificationQuality(Log10_Function_Output, Annotation)
NF = size(FunctionArray,1);

nAnnotation = ~Annotation;
E_pos = mean(Log10_Function_Output(:, Annotation),2); % This is the expected output level for the positive cases in each circuit
E_neg = mean(Log10_Function_Output(:,nAnnotation),2); % This is the expected output level for the negative cases in each circuit
AverageMargin = E_pos - E_neg;

% Re-weighting Samples
Worst_Pos = min(Log10_Function_Output(:, Annotation),[],2);
Worst_Neg = max(Log10_Function_Output(:,nAnnotation),[],2);

Clasification_pwr = zeros(size(Log10_Function_Output));
Clasification_pwr(:, Annotation) =  Log10_Function_Output(:, Annotation) - repmat(Worst_Neg, 1, sum( Annotation));
Clasification_pwr(:,nAnnotation) = -Log10_Function_Output(:,nAnnotation) + repmat(Worst_Pos, 1, sum(nAnnotation));
Clasification_pwr = max(Clasification_pwr,0);
MisClasifications = Clasification_pwr > 0;

SampleWeights = sum(MisClasifications,1) ./ NF;
SampleWeights = SampleWeights / sum(SampleWeights);
WeightedPerf  = Clasification_pwr * SampleWeights';
end