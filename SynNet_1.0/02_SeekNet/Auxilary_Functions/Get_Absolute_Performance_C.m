% Out put is the expected difference between the positive and negative
% inputs. In the case of boolean functions this is equivalent to: TPR - FPR
% where TPR is True Pos Rate, and FPR is False Pos Rate
% This is also equal to "informedness measure, which is the TPR + SPC - 1
% where SPC is Specifity or True Negative Rate
% Pejman

function    [AbsolutePerformance, Function_Output, Continuous_Stats] = Get_Absolute_Performance_C(FunctionArray, Data, Annotation)
% Function_Output = log10(Evaluate_Function_C(FunctionArray, Data));
% E_pos = mean(Function_Output(:, Annotation),2); % This is the expected output level for the positive cases in each circuit
% E_neg = mean(Function_Output(:,~Annotation),2); % This is the expected output level for the negative cases in each circuit
% AbsolutePerformance     = (E_pos - E_neg);


%% New performance (AUC + margin)
N = size(FunctionArray,1);
nAnnotation = ~Annotation;
Function_Output = log10(Evaluate_Function_C(FunctionArray, Data));

E_pos = mean(Function_Output(:, Annotation),2); % This is the expected output level for the positive cases in each circuit
E_neg = mean(Function_Output(:,nAnnotation),2); % This is the expected output level for the negative cases in each circuit
Continuous_Stats.MarginA = (E_pos - E_neg);

m_Pos = min(Function_Output(:, Annotation),[],2);
M_Neg = max(Function_Output(:,nAnnotation),[],2);
Continuous_Stats.MarginW = (m_Pos - M_Neg);

Continuous_Stats.Margin = Combined_Margins(Continuous_Stats);

dA = double(Annotation);
for i = N:-1:1
%     [~, ~, ~, AUC(i,1)] = perfcurve2(Annotation, Function_Output(i,:),  true);
    AUC(i,1) =    fastAUC(dA, Function_Output(i,:),  1);
end
AUC = AUC*100;
if nargout > 2
    Continuous_Stats.AUC     = AUC;
end

Normalized_Continuous_Margin = Continuous_Stats.Margin * 1E-6; % This means I assume margin would be much less than 1E+6, which given that it's in log10 scale it's fine, if you're uncomfortable with that use the following two commented lines instead!
AbsolutePerformance     = AUC + Normalized_Continuous_Margin;
% [~,I] = sortrows([AUC, MarginA], [1 2]);
% AbsolutePerformance(I,1) = linspace(0,AUC*1E-4*MarginA,N);



end