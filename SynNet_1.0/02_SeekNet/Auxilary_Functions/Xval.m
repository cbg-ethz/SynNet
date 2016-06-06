function [CV_Stats, StatsLabels] = Xval(AnalysisReportFile, CV_folds)
global Consts

if nargin==0
    AnalysisReportFile = '../F03_Results/BreastCancer_Data_C/C2-IDC.mat';
    CV_folds =3;
end

load(AnalysisReportFile);
fprintf('Cross-validating "%s"\n', Sim.CancerLbl)
Consts   = Sim.Consts;

%Consts.Learning_Convergence_Thr_1 = 100;
%Consts.Learning_Convergence_Thr_2 = 10;
% Consts.Plotting = false;
% Consts.Learning_MaxHitsToReport   = 1;

CV_JointReportName = [Consts.JointReportName(1:end-4) '-CV.txt'];

%% Cross-Validation
S0 = find(~Sim.Data.Annots);
S1 = find( Sim.Data.Annots);
L0 = length(S0);
L1 = length(S1);
mL = min(L0, L1);
if mL< CV_folds
    fprintf('%d-fold Cross-validation was cancelled due to lack of samples (one of the classes has only %d samples).\n',CV_folds, mL);
    CV_Stats = nan;
    StatsLabels = [];
    return
end

%% Make Cross-validation sets
B0 = floor(L0/CV_folds);
B1 = floor(L1/CV_folds);

PI0 = randperm(L0);
PI1 = randperm(L1);
% Break it to equal sets
for k = CV_folds:-1:1
    t0 = (B0*(k-1)+1):(B0*k);
    t1 = (B1*(k-1)+1):(B1*k);
    CV_S0{k} = S0(PI0(t0));
    CV_S1{k} = S1(PI1(t1));
end

% add the leftover samples in the sets till data is finished
for k = 1:CV_folds
    t = B0*CV_folds+k;
    if t>L0
        % its finished
        break;
    end
    CV_S0{k} = [CV_S0{k}, S0(PI0(t))];
end

% add the leftover samples in the sets till data is finished
for k = 1:CV_folds
    t = B1*CV_folds+k;
    if t>L1
        % its finished
        break;
    end
    CV_S1{k} = [CV_S1{k}, S1(PI1(t))];
end


CV_Stats = nan(CV_folds,6);
for k = CV_folds:-1:1
    t = [1:(k-1) k+1:CV_folds];
    CV_sim_L = SimSelect(Sim,[CV_S0{t} CV_S1{t}]);
    CV_sim_T = SimSelect(Sim,[CV_S0{k} CV_S1{k}]);
    
    %% Learn model with Learning set
    CV_sim_L = Fit_Classifier(CV_sim_L);
    
    %% Evaluate CV-Error of the Best function on the Test set
    Learned_Function = CV_sim_L.Results.BestCircuits(1,:,:);
    [CV_Stats(k,1:3), StatsLabels] = Get_AllStats(Learned_Function, CV_sim_T.Data);
    
    %% Evaluate CV-Error of the Pruned Best function on the Test set
    Learned_Function_Pruned = Prune_Circuit(Learned_Function, CV_sim_L.Data);
    [CV_Stats(k,4:6), ~] = Get_AllStats(Learned_Function_Pruned, CV_sim_T.Data);
    
    CV.TData{k} = CV_sim_T.Data;
    CV.Learned_Function{k} = Learned_Function;
    CV.Learned_Function_Pruned{k} = Learned_Function_Pruned;
end
%% Evaluate Original Error of the Best function on the Test set
Full_Data_Function = Sim.Results.BestCircuits(1,:,:);
[FL_Stats(1,1:3), ~] = Get_AllStats(Sim.Results.BestCircuits(1,:,:), Sim.Data);

Full_Data_Function_Pruned = Prune_Circuit(Full_Data_Function, Sim.Data);
[FL_Stats(1,4:6), ~] = Get_AllStats(Full_Data_Function_Pruned, Sim.Data);
StatsLabels = [StatsLabels strcat(StatsLabels, '_Prun')];

%% Make the report
% Check if Joiont report exists
JFisEmpty =  ~exist(CV_JointReportName, 'file');
Statsout = mean(CV_Stats,1);
StatsSEM = std(CV_Stats,1,1)/sqrt(size(CV_Stats,1)); % standard error of the mean
statsLout = [StatsLabels; strcat('CV_', [StatsLabels; strcat(StatsLabels, '_SE')])];

if JFisEmpty
    JFout = fopen(CV_JointReportName, 'w');
    fprintf(JFout, 'CancerType\tUsed_File');
    fprintf(JFout, '\t%s', statsLout{:});
    fprintf(JFout, '\tCV_folds\n');
else
    JFout = fopen(CV_JointReportName, 'a');
end

% Write the joint dataset summary report
fprintf(JFout, '%s\t%s', Sim.CancerLbl, FoutName(1:end-4));
fprintf(JFout, '\t%.2f', [FL_Stats; Statsout; StatsSEM]);
fprintf(JFout, '\t%d\n', CV_folds);
fclose(JFout);

% Add CV stats to the original file
load(AnalysisReportFile);
Sim.Results.CV       = CV;
Sim.Results.CV.Stats = CV_Stats;
Sim.Results.CV.StatsLabels = StatsLabels;
save(AnalysisReportFile, 'Sim', 'FoutName')
fprintf('Cross-validation done!\n')
end


function SubProblem = SimSelect(Problem, I)
% This just subselects some data points.
SubProblem = Problem;

D           = Problem.Data;
D.Values    = D.Values(:,I);
D.BValues   = D.BValues(:,I);
D.BMargin   = D.BMargin(:,I);

D.Labels    = D.Labels(I);
D.Annots    = D.Annots(I);
D.SampleID  = D.SampleID(I);


SubProblem.CancerLbl = Problem.CancerLbl;
SubProblem.Consts    = Problem.Consts;
SubProblem.Data      = D;
end