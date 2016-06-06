function MakeSimulationDataReport(FoutName, Sim)
global Consts
BestC = Sim.Results.BestCircuits(1,:,:);
SyntC = Sim.Synthesis_Function;
BS_int = intersect(BestC(:), SyntC(:)); BS_int(BS_int==0)=[];
[Sim.TestData.BValues, ~ , Sim.TestData.BMargin, Problem.TestData.Threshold] = Quantize_Expresison(Sim.TestData.Values);

% [~, Continous_MarginW, Continous_MarginA, Continous_AUC] = Get_Continous_Margins(BestC, Sim.Data.Values, Sim.Data.Annots);
% [~, Sim.Result.Output_C              ] = Get_Absolute_Performance_C(BestC, Sim.Data.Values, Sim.Data.Annots);
% [~, Sim.Result.Output_B, Binary_Stats] = Get_Absolute_Performance_B(BestC, Sim.Data       , Sim.Data.Annots);
[LernDstats, StatsLabels] = Get_AllStats(BestC, Sim.Data);
TestDstats                = Get_AllStats(BestC, Sim.TestData);
LernDstats_Synt           = Get_AllStats(SyntC, Sim.Data);
TestDstats_Synt           = Get_AllStats(SyntC, Sim.TestData);
if strcmpi(Consts.Learning_Pruning, 'off')
    % Skip pruning
    LernDstats_Prun           = nan(size(LernDstats));
    TestDstats_Prun           = nan(size(TestDstats));
    Out_tPrunCircuitSize      = nan;
    PS_int = [];
else
    % Prune the results
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

% [~, Continous_MarginW_t, Continous_MarginA_t, Continous_AUC_t] = Get_Continous_Margins(BestC, Sim.TestData.Values, Sim.TestData.Annots);
% [~, Sim.Result.Output_C_t              ]   = Get_Absolute_Performance_C(BestC, Sim.TestData.Values, Sim.TestData.Annots);
% [~, Sim.Result.Output_B_t, Binary_Stats_t] = Get_Absolute_Performance_B(BestC, Sim.TestData       , Sim.TestData.Annots);
%% Add a header line, if report is empty
Header = sprintf('CancerType\tGenes\tGenesinSynthCirc\tSamples\tPositiveSamples\tPositiveRate(%%)\tAnalysisMode\tHitMaxRound\tWasStuck\tGenesinLearntCirc\tSyn_and_Best\tGenesinPrunedCirc\tSyn_and_Prun\tTime(sec)');
Header = [Header,  sprintf('\t%s_L'   , StatsLabels{:})];
Header = [Header,  sprintf('\t%s_T'   , StatsLabels{:})];
Header = [Header,  sprintf('\t%s_Lsyn', StatsLabels{:})];
Header = [Header,  sprintf('\t%s_Tsyn', StatsLabels{:})];
Header = [Header,  sprintf('\t%s_Lpru', StatsLabels{:})];
Header = [Header,  sprintf('\t%s_Tpru', StatsLabels{:})];

if Sim.SimulationNumber==1
    JFout = fopen(FoutName, 'a');
    fprintf(JFout, '%s', Header);
    fprintf(JFout, '\n');
    fclose(JFout);
end


%% Preparing the output
Out_CancerType      = Sim.CancerLbl; % Name of the cancer type being analyzed.
Out_Genes           = size(Sim.Data.Values,1); % Number of genes included in the analysis.
Out_tSynCircuitSize = sum(SyntC(:)~=0); % number of genes in the circuit used for data synthesis
Out_Samples         = size(Sim.Data.Values,2); % Number of samples included in the analysis.
Out_PositiveSamples = sum(Sim.Data.Annots==true); % Number of positives samples in data
Out_PositiveRate    = Out_PositiveSamples/Out_Samples * 100; % Percentage of positives samples in data
Out_AnalysisMode    = Consts.AnalysisMode; % The used mode of analysis: "C"(Continuous) or "B"(Binary).
Out_HitMaxRound     = Sim.Results.ExitFlags.HitMaxRounds; % If the optimization stopped due to exhausting "Learning_MaxOptimizationRounds" constant.
Out_WasStuck        = Sim.Results.ExitFlags.WasStuck; % What was the longest strech of rounds that optimizer failed to improve the best classifier.
Out_Time            = Sim.Results.ExitFlags.ElapsedTime; % How long did the optimization take (in seconds)
Out_tBestCircuitSize= sum(BestC(:)~=0); % number of genes in the circuit used for data synthesis
Out_Syn_and_Best    = length(BS_int); % Number of genes shared between best and Synthetic network
Out_Syn_and_Prun    = length(PS_int); % Number of genes shared between best and Pruned network

%% Write Output file
RowBuff = [
    sprintf('%s\t', Out_CancerType),...
    sprintf('%d\t', Out_Genes, Out_tSynCircuitSize, Out_Samples, Out_PositiveSamples),...
    sprintf('%.2f\t', Out_PositiveRate),...
    sprintf('%c\t', Out_AnalysisMode),...
    sprintf('%d\t%d\t%d\t%d\t%d\t%d', Out_HitMaxRound, Out_WasStuck, Out_tBestCircuitSize, Out_Syn_and_Best, Out_tPrunCircuitSize, Out_Syn_and_Prun),...
    sprintf('\t%.2f', [Out_Time, LernDstats, TestDstats, LernDstats_Synt, TestDstats_Synt, LernDstats_Prun, TestDstats_Prun]),...
    sprintf('\n')];

% Write the joint dataset summary report
JFout = fopen(FoutName, 'a');
fprintf(JFout, '%s', RowBuff);
fclose(JFout);
end