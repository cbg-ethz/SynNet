% This file gets a gene expression matrix and an annotation and builds up a
% list including the genes that can annotate each specific sample
% correctly. These genes are then weighted by how good they are overall. 
% The script makes (2*number of genes) classifiers, and index them such
% that the first half correspond to single genes, and the second half correspond to the negated genes.

% Written by Pejman
% Jan 18th, 2013, Basel
% ------------------------------------------

function [PrimeNets] = BuildRecruitIndex_B(Problem)

BinData = Problem.Data.BValues; 
Annoation = Problem.Data.Annots; % Sample annotation, Cancer=1, healthy=0

Ngene     = size(BinData,1);

%% Pre-evaluate single gene classifiers
% Here I evaluate all signle gene classifiers, and build an index list,
% "RecruitList{s}", which says which genes can be used to correctly classify a
% specific sample "s".
NFhits = 2 * Ngene;
PrimeNets.Function  = Construct_Function([1:Ngene -(1:Ngene)]);
PrimeNets.IsaNOT    = [false(Ngene,1); true(Ngene,1)]; % Just denoting that this is a "NOT" circuit
PrimeNets.GeneID    = [(1:Ngene)'; (1:Ngene)'];

[AtNOutput, AtNOutputIsVague] = Evaluate_Function_B(PrimeNets.Function, BinData);
AtNet_perf = (AtNOutput == repmat(Annoation,NFhits, 1)); % Find misclassified outputs.
AtNet_perf(AtNOutputIsVague)  = false;       % Additionally, consider vague outputs as misclassified.

Sample_easiness = sum(AtNet_perf,1);
if any(Sample_easiness==0)
    warning(sprintf('Theres one or more unclassifiable sample: %s', SLabels{Sample_easiness==0}));
    Sample_easiness(Sample_easiness==0)= 1;
end
AtNet_perf = AtNet_perf./repmat(Sample_easiness, size(AtNet_perf,1),1); % Normalize the total sum to one.
OveralPerf = sum(AtNet_perf,2); % Overall performance of a gene.

% Sort the list so that most probable hits would be first, this will
% speed up the finding in the sampling of multinomial.
[PrimeNets.AtNet_Overalperf, PerfOrder] = sort(OveralPerf, 'descend');
PrimeNets.Function = PrimeNets.Function(PerfOrder,:,:);
PrimeNets.IsaNOT   = PrimeNets.IsaNOT(PerfOrder);
PrimeNets.GeneID   = PrimeNets.Function(:,1,1);  % This has negative values for the 
AtNet_perf = AtNet_perf(PerfOrder, :);
PrimeNets.SampleWise_WeightCumSum = cumsum(AtNet_perf,1);
PrimeNets.SampleWise_WeightCumSum(end,:) = 1; % Just to make sure total sum is exactly one, otherwise it might be tiny bit different due to numeric errors.
end





