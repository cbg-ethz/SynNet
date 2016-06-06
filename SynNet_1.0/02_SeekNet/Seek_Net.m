function Seek_Net(Constraints)
global Consts
%% Initialization
set(0, 'defaultTextInterpreter', 'tex');
try addpath('Auxilary_Functions/');     catch ; end

if nargin <1; Constraints = 'Constraints.txt'; end
Consts = Fetch_Constraints(Constraints);


Cnt= 0;
if Consts.IsSimulation
    %% Simulation Analysis
    FoutName = Initialize_Simulation(Constraints);
    while Cnt< Consts.Simulation_Rounds
        Cnt = Cnt+1;
        fprintf('\n[Simulation Round %d] - ', Cnt);
        try
            Sim = [];
            while isempty(Sim)
                Sim = Simulate_data();
            end
            Sim.CancerLbl = ['Sim-C' int2str(Cnt)];
            
            Sim = Fit_Classifier(Sim);
            Sim.SimulationNumber = Cnt;
            MakeSimulationDataReport(FoutName, Sim);
            
        catch Err
            disp('Something went wrong by the way!')
            disp(Err.getReport)
        end
        fprintf('\n\n')
    end
else
    %% Real data
    mkdir(Consts.OutFolder);
    [miRNA_Data, MetaData] = Read_miRNAData(Consts.DataPath, Consts.MetaDataPath);
    fclose(fopen(Consts.JointReportName, 'w')); % Make an empty file for the joint summary report
    Cid = unique(miRNA_Data.Annotation); Cid(Cid<1)=[];
    
    for tmpCancerIndex = 1:length(Cid)%:-1:1
        CancerIndex = Cid(tmpCancerIndex);
        CancerLbl = miRNA_Data.SampleLabels{find(miRNA_Data.Annotation == CancerIndex, 1, 'first')};
        fprintf('\n[Analysing - %s]\n', CancerLbl);
        FoutName = [Consts.OutFolder '/' Consts.AnalysisMode int2str(CancerIndex) '-' CancerLbl '.txt'];
        copyfile(Constraints,[FoutName(1:end-4) '-Used_Constants.log'],'f');
        
        %         try
        Sim = miRNA_Feat_Filt(miRNA_Data, 0, CancerIndex);
        Sim.CancerLbl = CancerLbl;
        Sim = Fit_Classifier(Sim);
        Sim.MetaData = MetaData;
        if strcmpi(Consts.Learning_Pruning, 'off')
            % Make only one report
            MakeRealDataReport(FoutName, Sim);
        else
            % Make one extra report with pruned results
            Consts.Learning_Pruning = 'off';
            MakeRealDataReport(FoutName, Sim);
            Consts.Learning_Pruning = 'on';
            MakeRealDataReport([FoutName(1:end-4) '-Prune.txt'], Sim);
            % Do 3-fold X-validation for Normal and Pruned Networks in on go
            if strcmpi(Consts.Learning_CrossValidation, 'on')
                Xval([FoutName(1:end-4) '.mat'], 3);
            end
        end
        %         catch Err
        %             disp('Something went wrong by the way!'); disp(Err.getReport)
        %         end
        fprintf('\n\n')
    end
end
end

