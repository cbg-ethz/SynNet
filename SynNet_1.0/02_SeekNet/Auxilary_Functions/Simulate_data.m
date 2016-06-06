% This code simulates data.
% output is an structure to be used for the Seek_DoomNet.m
% This file is called from the "Seek_DoomNet.m" if it's in the synthetic
% data mode.
%
% Written By Pejman, 15Jan2013, Basel
% Pejman.m@gmail.com
%----------------------------------------------

function Sim = Simulate_data()
global Consts
addpath('Auxilary_Functions/')
Synthesis_ratio = Consts.Simulation_DataSynthesis_ratio;% % If this ratio is increased from 1
Mu       = log10(Consts.Quantizer_Threshold);
Ngene    = randi(Consts.Simulation_Max_Genes  - Consts.Simulation_Min_Genes   +1) + Consts.Simulation_Min_Genes   -1;
Nsamples = randi(Consts.Simulation_Max_Samples- Consts.Simulation_Min_Samples +1) + Consts.Simulation_Min_Samples -1;

Nsynth  = Nsamples * 2 * Synthesis_ratio; % We make twice as much samples because later we make one leanset and one testset.
% Function is encoded such that
% A number shows a gene (literals).
% Each row corresponds to an OR gate
% Column are ANDed together after the row-wise ORs.
% Negative numbers correspond to Negated literals.
Function = Build_Random_Function(Ngene, Consts);
% while length(Annoation)<Nsamples
%% Data simulation
% fix the simulation seed (For debugging only1)
% s = RandStream('mcg16807','Seed',0);
% try
% RandStream.setDefaultStream(s);
% catch tmpE
%     RandStream.setGlobalStream(s);
% end
Ldata = randn(Ngene, Nsynth);% produce log-data from standard normal samples.
ScalingFactor = rand(Ngene,1) - .5;%3 / quantile(abs(Ldata(:)), .95) ; % SO 95% of the point fall within 3 folds far from the step %(log(ZeroLvl) - Mu) / min(log(1-eps),quantile(Ldata(:), (1 - VagueFreq)/2));% scale data to tune the variance to get roughly the right portion of vague values
MuG = rand(Ngene, 1)*2 + Mu - 1;
tmpGenesData = 10.^(Ldata .* repmat(ScalingFactor, 1, Nsynth) + repmat(MuG, 1, Nsynth));

[BValues, ~, BMargin] = Quantize_Expresison(tmpGenesData);
[Function_Output, Margin] = Evaluate_Function_B_Margin(Function, BValues, BMargin);
MThr = quantile(Margin, 1- 1/Synthesis_ratio);
MFilt= Margin>=MThr;
if Consts.AnalysisMode == 'B'
    tmpAnnot = Function_Output;
    Exclude = ~MFilt;
else
    Npos = sum( Function_Output(MFilt));
    Nneg = sum(~Function_Output(MFilt));

    Function_Output_C = log10(Evaluate_Function_C(Function, tmpGenesData));
    [~, I] = sort(Function_Output_C, 'descend');
    
    tmpAnnot = nan(size(Function_Output_C));
    tmpAnnot(I(1:Npos)) =  1;
    tmpAnnot(I((end-Nneg+1):end)) =  0;
    Exclude = isnan(tmpAnnot);
end

for i = Ngene:-1:1
    GenIDs{i} = ['G' int2str(i)];
end

% Remove Vague samples
tmpAnnot(Exclude) = [];
tmpGenesData(:,Exclude) = [];
% Function_Output_C(Exclude) = [];

Annoation = [ tmpAnnot == 1];
GenesData = [ tmpGenesData];

for i = length(Annoation):-1:1
    if Annoation(i)==true
        Labels{i} = ['Cancer' int2str(i)];
    else
        Labels{i} = ['Healthy' int2str(i)];
    end
end


allData.Values  = GenesData;
allData.GenIDs  = GenIDs;
allData.Labels  = Labels;
allData.Annots  = Annoation; % Zeros for healthy, ones for tumor
allData.SampleID= 1:(Nsamples*2);

TestSamples = randperm(Nsamples*2, Nsamples); % we choose half of the data as test set

LearnD = allData;
LearnD.Values(:,TestSamples)=[];
LearnD.Annots(TestSamples)  =[];
if length(unique(LearnD.Annots))== 1
    % All samples are only positive, or only negative, discard the whole
    % thing
    Sim = [];
    return
end
Sim.Data = LearnD; clear LearnD

TestD = allData; clear allData
TestD.Values(:,~TestSamples)=[];
TestD.Annots(~TestSamples)  =[];
Sim.TestData = TestD; clear TestD

Sim.Synthesis_Function = Function;

fprintf('Synthetic data: %d genes, %d samples, and %d genes in target circuit.\n', Ngene, Nsamples, sum(abs(Function(:))>0))
fprintf('================================\n')
end

