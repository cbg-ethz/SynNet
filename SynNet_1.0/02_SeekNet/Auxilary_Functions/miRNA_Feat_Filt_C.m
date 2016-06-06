% This code does couple ofstuff:
%       1: Quantizes the data to Zero, One, and NaN(vague, i.e. between zero and one)  based on ZeroLvl, and
%       OneLvl thresholds. 
%       2: Remove uninformative genes; those who are constant along all samples.
%       3: Reoprts if the are some unresolvable samples; those that have all genes(after excluding non-informative genes) in NaN state. 
%

% Pejman Jun. 2014, Pejman.m@gmail.com
%----------------------------------

function Sim = miRNA_Feat_Filt_C(miRNA_Data, HealthyIndex, TumorIndex)
global Consts

%% Normalize Gene Expressions
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

if isfield(Consts, 'MiRNAs_Blacklist') && ~isempty(Consts.MiRNAs_Blacklist)
    % Exclude these miRNas from the analysis
    miRNA_Data = Discard_Blacklisted_miRNAs(miRNA_Data, Consts);
end




miRNA_Data.Dnormed = miRNA_Data.Dnormed + miRNA_Data.Pseudocount;% We add pseudo-count, here at the very end, after normalizing for the library size 

%% Remove genes with too little effect
ExpRanges(:,1) = min(miRNA_Data.Dnormed(:,All),[],2);
ExpRanges(:,2) = max(miRNA_Data.Dnormed(:,All),[],2);
BaseFuncs = zeros(2, Consts.MaxAnd, Consts.MaxOr);
BaseFuncs(1,1,1) = +1;  % A basic literal
BaseFuncs(2,1,1) = -1;  % A basic negated literal

for i = size(ExpRanges,1):-1:1
    tmpReponse = log2(Evaluate_Function_C(BaseFuncs, ExpRanges(i,:)));
    ResponseRange(i,1) = max(abs(tmpReponse(:,2)-tmpReponse(:,1))); % this is the most effect a gene can make alone
end

UnInfFilt = ResponseRange < Consts.Learning_mRT_log2 ;
fprintf('%d out of %d genes were excluded from analysis due to "Learning_mRT_log2" threshold\n', sum(UnInfFilt), length(UnInfFilt));


%% Make output
Sim.Data.Values  = miRNA_Data.Dnormed(~UnInfFilt,All);
Sim.Data.GenIDs  = miRNA_Data.GeneNames(~UnInfFilt);
Sim.Data.Labels  = miRNA_Data.SampleLabels(All);
Sim.Data.Annots  = AllAnnot; % Zeros for healthy, ones for tumor
Sim.Data.SampleID= miRNA_Data.SampleID(All);
Sim.Data.GeneGroupNames = miRNA_Data.GeneGroupNames(~UnInfFilt);
Sim.Data.Mature_Sequence = miRNA_Data.Mature_Sequence(~UnInfFilt);
end