% This code does couple ofstuff:
%       1: Quantizes the data to Zero, One, and NaN(vague, i.e. between zero and one)  based on ZeroLvl, and
%       OneLvl thresholds.
%       2: Remove uninformative genes; those who are constant along all samples.
%       3: Reoprts if the are some unresolvable samples; those that have all genes(after excluding non-informative genes) in NaN state.
%

% Pejman Jun. 2014, Pejman.m@gmail.com
%----------------------------------

function Sim = miRNA_Feat_Filt_B(miRNA_Data, HealthyIndex, TumorIndex)
global Consts

%% Uninformative Feat Filtering
ZeroLvl = Consts.Quantizer_ZeroLvl;
OneLvl  = Consts.Quantizer_OneLvl;

Healthy = find(miRNA_Data.Annotation==HealthyIndex);
Tumors  = find(miRNA_Data.Annotation==TumorIndex);

All     = [Healthy' Tumors'];
AllAnnot= [false(1,length(Healthy)) true(1, length(Tumors))]; % Zeros for healthy, ones for tumor

miRNA_Data = normalize_miRNAs(miRNA_Data, All, Consts);
if isfield(Consts, 'Mature_miRNAs') && ~isempty(Consts.Mature_miRNAs)
    % Use the given miRNA sequences to combine the expressions of similar miRNAs
    miRNA_Data = Combine_Similar_miRNAs(miRNA_Data, Consts);
else
    miRNA_Data.GeneGroupNames = cell(size(miRNA_Data.GeneNames));
    miRNA_Data.Mature_Sequence = cell(size(miRNA_Data.GeneNames));
end
miRNA_Data.Dnormed = miRNA_Data.Dnormed + miRNA_Data.Pseudocount;% We add one pseudo-count, here at the very end, after normalizing for the library size

D = miRNA_Data.Dnormed;

Dbin = nan(size(D));

Dbin(D<=ZeroLvl) = 0;
Dbin(D>=OneLvl)  = 1;
%[BinData, BinDataIsVague] = Quantize_Expresison(GenesData, ZeroLvl, OneLvl)

UnInfFilt = ...
    sum(isnan(Dbin(:,All)   )| Dbin(:,All)==0 ,2) == length(All)     | ...
    sum(isnan(Dbin(:,All)   )| Dbin(:,All)==1 ,2) == length(All)     ;

UnInfSamples = ...
    sum(isnan(Dbin(~UnInfFilt,All)),1) == sum(~UnInfFilt)     ;

if sum(UnInfSamples) > 0
    beep
    disp('WARNING: There are unresolvable samples in the dataset:');
    disp(miRNA_Data.SampleLabels(All(UnInfSamples)));
end



Sim.Data.Values  = miRNA_Data.Dnormed(~UnInfFilt,(All(~UnInfSamples)));
Sim.Data.GenIDs  = miRNA_Data.GeneNames(~UnInfFilt);
Sim.Data.Labels  = miRNA_Data.SampleLabels(All(~UnInfSamples));
Sim.Data.Annots  = AllAnnot(~UnInfSamples); % Zeros for healthy, ones for tumor
Sim.Data.SampleID= miRNA_Data.SampleID(All(~UnInfSamples));
Sim.Data.GeneGroupNames = miRNA_Data.GeneGroupNames(~UnInfFilt);
Sim.Data.Mature_Sequence = miRNA_Data.Mature_Sequence(~UnInfFilt);
end