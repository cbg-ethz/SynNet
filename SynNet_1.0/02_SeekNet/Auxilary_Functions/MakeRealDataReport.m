function MakeRealDataReport(FoutName, Sim)
global Consts
if nargin == 1
    % then Assume this is a '.mat' file that has both Sim and FoutName in
    % it.
    load(FoutName, 'Sim');
    Consts   = Sim.Consts;
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
    Sim.Results.BestCircuits   = Prune_Circuit(Sim.Results.BestCircuits, Sim.Data);
end

[~, Continous_MarginW, Continous_MarginA, Continous_AUC] = Get_Continous_Margins(Sim.Results.BestCircuits(1,:,:), Sim.Data.Values, Sim.Data.Annots, [FoutName(1:end-4) '-OutLvl_C']);
[~, Sim.Result.Output_C, Contin_Stats] = Get_Absolute_Performance_C(Sim.Results.BestCircuits, Sim.Data.Values, Sim.Data.Annots);
[~, Sim.Result.Output_B, Binary_Stats] = Get_Absolute_Performance_B(Sim.Results.BestCircuits, Sim.Data       , Sim.Data.Annots);
[~, Sim.Result.Output_B_Margins] = Evaluate_Function_B_Margin(Sim.Results.BestCircuits, Sim.Data.BValues, Sim.Data.BMargin);
if Consts.AnalysisMode == 'B'
    Write_Function_InOuts_B(Sim.Results.BestCircuits(1,:,:), Sim.Data, Pej_Struct_RowSelect(Sim.Result, 1), [FoutName(1:end-4) '-Inputs']);
    Plot_binary_output(Sim.Results.BestCircuits(1,:,:), Sim.Data, Sim.Data.Annots, [FoutName(1:end-4) '-OutLvl_B'])
else
    Write_Function_InOuts_C(Sim.Results.BestCircuits(1,:,:), Sim.Data, Pej_Struct_RowSelect(Sim.Result, 1), [FoutName(1:end-4) '-Inputs']);
end
%% Preparing the output
Out_CancerType      = Sim.CancerLbl; % Name of the cancer type being analyzed.
Out_Genes           = size(Sim.Data.Values,1); % Number of genes included in the analysis.
Out_Samples         = size(Sim.Data.Values,2); % Number of samples included in the analysis.
Out_PositiveSamples = sum(Sim.Data.Annots==true); % Number of positives samples in data
Out_PositiveRate    = Out_PositiveSamples/Out_Samples * 100; % Percentage of positives samples in data
Out_AnalysisMode    = Consts.AnalysisMode; % The used mode of analysis: "C"(Continuous) or "B"(Binary).
Out_Pruning         = Consts.Learning_Pruning;
Out_HitMaxRound     = Sim.Results.ExitFlags.HitMaxRounds; % If the optimization stopped due to exhausting "Learning_MaxOptimizationRounds" constant.
Out_WasStuck        = Sim.Results.ExitFlags.WasStuck; % What was the longest strech of rounds that optimizer failed to improve the best classifier.
Out_Time            = Sim.Results.ExitFlags.ElapsedTime; % How long did the optimization take (in seconds)
Out_C_MarginA       = Continous_MarginA; % Average continuous classification margin of the best circuit (in log10 scale).
Out_C_MarginW       = Continous_MarginW; % Worst   continuous classification margin of the best circuit (in log10 scale).
Out_C_AUC           = Continous_AUC; % This is the Area under RoC curve of the continuous classification. It says how good the overall ranking is.
Out_B_Performed     = Binary_Stats.RawPerformance(1)*100; % Classification AUC of binary classification using the best circuit.
Out_B_MarginA       = Binary_Stats.MarginA(1);
Out_B_MarginW       = Binary_Stats.MarginW(1);
Out_B_FPrate        = Binary_Stats.FPrate(1); % False Positive rate of binary classification using the best circuit.
Out_B_FNrate        = Binary_Stats.FNrate(1); % False Negative rate of binary classification using the best circuit.

%% Write Output file
Header = sprintf('AnalysisName\tGenes\tSamples\tPositiveSamples\tPositiveRate(%%)\tAnalysisMode\tPruning\tHitMaxRound\tWasStuck\tTime(sec)\tC_MarginA\tC_MarginW\tC_AUC\tB_MarginA\tB_MarginW\tB_AUC\tB_FPrate\tB_FNrate\n');
RowBuff = [
    sprintf('%s\t', Out_CancerType),...
    sprintf('%d\t', Out_Genes, Out_Samples, Out_PositiveSamples),...
    sprintf('%.2f\t', Out_PositiveRate),...
    sprintf('%c\t', Out_AnalysisMode),...
    sprintf('%s\t', Out_Pruning),...
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
    BestsPerfs = [Binary_Stats.RawPerformance*100, Binary_Stats.MarginA];
else
    BestsPerfs = [Contin_Stats.AUC, Contin_Stats.MarginA];
end

Fout = fopen(FoutName, 'w');
fprintf(Fout, '%s', Header, RowBuff);
fprintf(Fout,'\nBest Circuit Found:\n%s\n', SPrint_Function(Sim.Results.BestCircuits, Sim.Data.GenIDs, BestsPerfs));


clear OutBuff
OutBuff{1}= [sprintf('Tag\tDetail1\tDetail2\tDetail3')];
OutBuff{end+1} = [sprintf('Sample Info\tSample IDs\t-\t-')     sprintf('\t%s', Sim.Data.SampleID{:})];
OutBuff{end+1} = [sprintf('Sample Info\tSample Names\t-\t-')   sprintf('\t%s', Sim.Data.Labels{:})];
OutBuff{end+1} = [sprintf('Sample Info\tDesired Annotation\t-\t(Binary)')   sprintf('\t%d', double(Sim.Data.Annots))];

for i = 1:Sim.Consts.Learning_MaxHitsToReport
    if Consts.AnalysisMode == 'B'
        OutBuff{end+1} = [sprintf('Circuit Info\tCircuit Output\tCircuit%d\t(Binary)', i)   sprintf('\t%d', double(Sim.Result.Output_B(i,:)))];
        OutBuff{end+1} = [sprintf('Circuit Info\tCircuit Output Margin\tCircuit%d\t(Ratio)', i)   sprintf('\t%.2f', 10.^Sim.Result.Output_B_Margins(i,:))];
    else
        OutBuff{end+1} = [sprintf('Circuit Info\tCircuit Output\tCircuit%d\t(mols/cell)', i)   sprintf('\t%.2f', 10.^Sim.Result.Output_C(i,:))];
    end
end
OutBuff{end+1} = sprintf('');
try
    % add meta data if available
    MD = Sim.MetaData.Values;
    iID = find(strcmp(Sim.MetaData.Labels, 'UniqueID'));
    [~, ai, bi]= intersect(Sim.Data.SampleID, Sim.MetaData.Values(:,iID));
    MD(ai,:)=Sim.MetaData.Values(bi,:);
    for i = 1:length(Sim.MetaData.Labels)
        if i == iID; continue;end
        OutBuff{end+1} = [sprintf('Sample Info\tSample Metadata\t%s\t-', Sim.MetaData.Labels{i}) sprintf('\t%s', MD{:,i})];
    end
catch
end

OutBuff{end+1} = sprintf('');
% add gene data
% OutBuff{end+1}= [sprintf('Tag\tmiRNA Name\tmiRNA Sequence\tmiRNA Group') ExpEmpty];

% Gene used in one of the reported circuit
CircuitGeneIDs = unique(abs(Sim.Results.BestCircuits(:)));% genes in all reported circiuts
CircuitGeneIDs(CircuitGeneIDs==0)=[];
GeneOB = [Sim.Data.GenIDs(CircuitGeneIDs) num2cell(Sim.Data.Values(CircuitGeneIDs,:))]';
CLOB = length(OutBuff);
OutBuff{CLOB+length(CircuitGeneIDs)}=[];
for i = 1:length(CircuitGeneIDs)
    j = CircuitGeneIDs(i);
    if isempty(Sim.Data.GeneGroupNames{j})
        RepSeq_= 'NA';
        RepGrn_= 'NA';
    else
        RepSeq_= Sim.Data.Mature_Sequence{j};
        RepGrn_= Sim.Data.GeneGroupNames{j};
    end
    OutBuff{CLOB+i}= [sprintf('In-Circuit miRNAs\t%s\t%s\t%s',...
        Sim.Data.GenIDs{j},...
        RepSeq_,...
        RepGrn_),...
        sprintf('\t%.2f', Sim.Data.Values(j,:)) ];
end

% Other genes included in the analysis
OtherGenes= 1:length(Sim.Data.GenIDs);
OtherGenes(CircuitGeneIDs)=[];
CLOB = length(OutBuff);
for i = 1:length(OtherGenes)
    j = OtherGenes(i);
    if isempty(Sim.Data.GeneGroupNames{j})
        RepSeq_= 'NA';
        RepGrn_= 'NA';
    else
        RepSeq_= Sim.Data.Mature_Sequence{j};
        RepGrn_= Sim.Data.GeneGroupNames{j};
    end
    OutBuff{CLOB+i}= [sprintf('Other miRNAs\t%s\t%s\t%s',...
        Sim.Data.GenIDs{j},...
        RepSeq_,...
        RepGrn_),...
        sprintf('\t%.2f', Sim.Data.Values(j,:)) ];
end

fprintf(Fout, '%s\n', OutBuff{:});
fclose(Fout);
end