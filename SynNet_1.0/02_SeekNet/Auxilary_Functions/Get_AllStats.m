function [Stats, StatsLabels] = Get_AllStats(Function, Data)
global Consts
% % This is how the stats are ordered
% Out_C_MarginA       = Continous_MarginA; % Average continuous classification margin of the best circuit (in log10 scale).
% Out_C_MarginW       = Continous_MarginW; % Worst   continuous classification margin of the best circuit (in log10 scale).
% Out_C_AUC           = Continous_AUC; % This is the Area under RoC curve of the continuous classification. It says how good the overall ranking is.
% Out_B_Performed     = Binary_Stats.RawPerformance*100; % Classification informedness of binary classification using the best circuit.
% Out_B_MarginA       = Binary_Stats.MarginA;
% Out_B_MarginW       = Binary_Stats.MarginW;

if Consts.AnalysisMode == 'B'
    
    
    [~, ~, Binary_Stats] = Get_Absolute_Performance_B(Function, Data , Data.Annots);
    Informedness  = Binary_Stats.RawPerformance; % Classification informedness of binary classification using the best circuit.
    Performed     = (Informedness +1) / 2; % This is equal to AUC, or to (Sensitivity + Specificity)/2
    MarginA       = Binary_Stats.MarginA;
    MarginW       = Binary_Stats.MarginW;
else
    [~, Continous_MarginW, Continous_MarginA, Continous_AUC] = Get_Continous_Margins(Function, Data.Values, Data.Annots);
    MarginA       = Continous_MarginA; % Average continuous classification margin of the best circuit (in log10 scale).
    MarginW       = Continous_MarginW; % Worst   continuous classification margin of the best circuit (in log10 scale).
    Performed     = Continous_AUC; % This is the Area under RoC curve of the continuous classification. It says how good the overall ranking is.
end
    Stats               = [MarginA, MarginW, 100 * Performed];
    StatsLabels  = {'MarginA' 'MarginW' 'AUC'};
end