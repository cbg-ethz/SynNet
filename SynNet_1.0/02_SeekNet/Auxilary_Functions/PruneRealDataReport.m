function PruneRealDataReport(FoutName, Sim)
global Consts
if nargin == 1
    % then Assume this is a '.mat' file that has both Sim and FoutName in
    % it.
    load(FoutName, 'Sim');
    Consts   = Sim.Constraints;
    FoutName = strrep(FoutName, '.mat', '');
    FoutName = [FoutName '.txt'];
    Consts.JointReportName = '';
else
    % Save the inputs in case they are needed later on
    save([FoutName(1:end-4) '.mat'], 'Sim', 'FoutName');
end


if strcmpi(Consts.Learning_Pruning, 'off')
    % Skip pruning
else
    % Prune the results
    for i 
    PrunC                     = Prune_Circuit(BestC, Sim.Data);
    Out_tPrunCircuitSize      = sum(PrunC(:)~=0);
    PS_int = intersect(PrunC(:), SyntC(:)); PS_int(PS_int==0)=[];

    if Out_tPrunCircuitSize == 0
        LernDstats_Prun           = [0 0 50];
        TestDstats_Prun           = [0 0 50];
    else
        LernDstats_Prun           = Get_AllStats(PrunC, Sim.Data);
        TestDstats_Prun           = Get_AllStats(PrunC, Sim.TestData);
    end
end

[~, Continous_MarginW, Continous_MarginA, Continous_AUC] = Get_Continous_Margins(Sim.Results.BestCircuits(1,:,:), Sim.Data.Values, Sim.Data.Annots, [FoutName(1:end-4) '-OutLvl_C']);
Plot_binary_output(Sim.Results.BestCircuits(1,:,:), Sim.Data, Sim.Data.Annots, [FoutName(1:end-4) '-OutLvl_B'])
[~, Sim.Result.Output_C              ] = Get_Absolute_Performance_C(Sim.Results.BestCircuits(1,:,:), Sim.Data.Values, Sim.Data.Annots);
[~, Sim.Result.Output_B, Binary_Stats] = Get_Absolute_Performance_B(Sim.Results.BestCircuits(1,:,:), Sim.Data       , Sim.Data.Annots);
Write_Function_InOuts(Sim.Results.BestCircuits(1,:,:), Sim.Data, Sim.Result, [FoutName(1:end-4) '-Inputs']);
%% Preparing the output
Out_CancerType      = Sim.CancerLbl; % Name of the cancer type being analyzed.
Out_Genes           = size(Sim.Data.Values,1); % Number of genes included in the analysis.
Out_Samples         = size(Sim.Data.Values,2); % Number of samples included in the analysis.
Out_PositiveSamples = sum(Sim.Data.Annots==true); % Number of positives samples in data
Out_PositiveRate    = Out_PositiveSamples/Out_Samples * 100; % Percentage of positives samples in data
Out_AnalysisMode    = Consts.AnalysisMode; % The used mode of analysis: "C"(Continuous) or "B"(Binary).
Out_HitMaxRound     = Sim.Results.ExitFlags.HitMaxRounds; % If the optimization stopped due to exhausting "Learning_MaxOptimizationRounds" constant.
Out_WasStuck        = Sim.Results.ExitFlags.WasStuck; % What was the longest strech of rounds that optimizer failed to improve the best classifier.
Out_Time            = Sim.Results.ExitFlags.ElapsedTime; % How long did the optimization take (in seconds)
Out_C_MarginA       = Continous_MarginA; % Average continuous classification margin of the best circuit (in log10 scale).
Out_C_MarginW       = Continous_MarginW; % Worst   continuous classification margin of the best circuit (in log10 scale).
Out_C_AUC           = Continous_AUC; % This is the Area under RoC curve of the continuous classification. It says how good the overall ranking is.
Out_B_Performed     = Binary_Stats.RawPerformance*100; % Classification AUC of binary classification using the best circuit.
Out_B_MarginA       = Binary_Stats.MarginA;
Out_B_MarginW       = Binary_Stats.MarginW;
Out_B_FPrate        = Binary_Stats.FPrate; % False Positive rate of binary classification using the best circuit.
Out_B_FNrate        = Binary_Stats.FNrate; % False Negative rate of binary classification using the best circuit.

%% Write Output file
Header = sprintf('CancerType\tGenes\tSamples\tPositiveSamples\tPositiveRate(%%)\tAnalysisMode\tHitMaxRound\tWasStuck\tTime(sec)\tC_MarginA\tC_MarginW\tC_AUC\tB_MarginA\tB_MarginW\tB_AUC\tB_FPrate\tB_FNrate\n');
RowBuff = [
    sprintf('%s\t', Out_CancerType),...
    sprintf('%d\t', Out_Genes, Out_Samples, Out_PositiveSamples),...
    sprintf('%.2f\t', Out_PositiveRate),...
    sprintf('%c\t', Out_AnalysisMode),...
    sprintf('%d\t', Out_HitMaxRound, Out_WasStuck),...
    sprintf('%.2f\t', Out_Time, Out_C_MarginA, Out_C_MarginW, Out_C_AUC, Out_B_MarginA, Out_B_MarginW, Out_B_Performed, Out_B_FPrate, Out_B_FNrate),...
    sprintf('\n')];

if ~isempty(Consts.JointReportName)
    % Check if Joiont report is empty
    JFout = fopen(Consts.JointReportName, 'r');
    fscanf(JFout,'%c',1);
    JFisEmpty =  feof(JFout);
    fclose(JFout);
    
    % Write the joint dataset summary report
    JFout = fopen(Consts.JointReportName, 'a');
    if JFisEmpty;    fprintf(JFout, '%s', Header); end
    fprintf(JFout, '%s', RowBuff);
    fclose(JFout);
end


if Consts.AnalysisMode == 'B'
    BestsPerfs = [Out_B_Performed, Out_B_MarginA];
else
    BestsPerfs = [Out_C_AUC, Out_C_MarginA];
end

Fout = fopen(FoutName, 'w');
fprintf(Fout, '%s', Header, RowBuff);
fprintf(Fout,'Best Circuit Found:\n\n%s\n', SPrint_Function(Sim.Results.BestCircuits, Sim.Data.GenIDs, BestsPerfs));
% fprintf(Fout, 'Samples failed to be classified:');
% fprintf(Fout,'\n%s', Sim.Results.BestCircuitsMissed{:});
%
% fprintf(Fout, '\n================================\nNumber of informative genes: %d\n', length(Sim.Data.GenIDs));
% fprintf(Fout, '\nQuantized miRNA Expression Composition\nSampleName\tHigh(%%)\tLow(%%)\tVague(%%)\n');
% Qstats = [sum(Sim.Results.Quantized_Data==1,1); sum(Sim.Results.Quantized_Data==0,1); sum(isnan(Sim.Results.Quantized_Data),1)]/length(Sim.Data.GenIDs) *100;
% OutBuff = [Sim.Data.Labels; num2cell(Qstats)];
% fprintf(Fout, '%s\t%.1f\t%.1f\t%.1f\n', OutBuff{:,:});
fclose(Fout);
% Pej_Write_Table([FoutName(1:end-4) '_QData.dat'], Sim.Results.Quantized_Data , Sim.Data.Labels, Sim.Data.GenIDs);
end