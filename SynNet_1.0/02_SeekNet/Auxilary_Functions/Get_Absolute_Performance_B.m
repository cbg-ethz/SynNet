% Out put is the expected difference between the positive and negative
% inputs. In the case of boolean functions this is equivalent to: TPR - FPR
% where TPR is True Pos Rate, and FPR is False Pos Rate
% This is also equal to "informedness measure, which is the TPR + SPC - 1
% where SPC is Specifity or True Negative Rate
% Pejman

function    [AbsolutePerformance, Function_Output, Binary_Stats] = Get_Absolute_Performance_B(FunctionArray, Data, Annotation)
global Consts

nAnnotation = ~Annotation;

[Function_Output, Output_Margins] = Evaluate_Function_B_Margin(FunctionArray, Data.BValues, Data.BMargin);
isNotCorrect = Function_Output ~= repmat(Annotation, size(Function_Output,1), 1);
Output_Margins(isNotCorrect) = -Output_Margins(isNotCorrect);
Binary_Margin_pos = mean(Output_Margins(:, Annotation),2);
Binary_Margin_neg = mean(Output_Margins(:,nAnnotation),2);
Binary_Stats.MarginA = (Binary_Margin_pos + Binary_Margin_neg)/2;
Binary_Stats.MarginW = min(Output_Margins,[],2);
Binary_Stats.Margin = Combined_Margins(Binary_Stats);

Normalized_Binary_Margin = Binary_Stats.Margin * Consts.BMargin_Coeficient; % To make sure it's always smaller that the smallest possible change due to change in classifier accuracy
if nargout > 2
    % The false positive rate is FP / (FP + TN).
    Binary_Stats.FPrate = sum(isNotCorrect(:,nAnnotation),2) / sum(nAnnotation)*100;
    % The false positive rate is FN / (FN + TP).
    Binary_Stats.FNrate = sum(isNotCorrect(:, Annotation),2) / sum( Annotation)*100;
end

E_pos = mean(Function_Output(:, Annotation),2); % This is the expected output level for the positive cases in each circuit
E_neg = mean(Function_Output(:,nAnnotation),2); % This is the expected output level for the negative cases in each circuit
AbsolutePerformance     = (E_pos - E_neg) + Normalized_Binary_Margin;

if nargout > 2
    Binary_Stats.RawPerformance = (E_pos - E_neg + 1)/2;% This is equal to AUC, or to (Sensitivity + Specificity)/2
end
end