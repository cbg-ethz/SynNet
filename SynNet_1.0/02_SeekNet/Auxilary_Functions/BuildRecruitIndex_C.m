% This file gets a gene expression matrix and an annotation and builds up a
% list based on how good a genes that can classify each specific sample.
% These genes are also weighted by how good they are overall. 
% The script makes (2*number of genes) classifiers, and index them such
% that the first half correspond to single genes, and the second half correspond to the negated genes.

% Written by Pejman
% Jan 18th, 2013, Basel
% ------------------------------------------

function [PrimeNets] = BuildRecruitIndex_C(Problem)
global Consts

GenesData = Problem.Data.Values; % Data matrix, with samples on rows, and genes on columns
Annotation = Problem.Data.Annots; % Sample annotation, Cancer=1, healthy=0
SLabels   = Problem.Data.Labels; % Sample Names
Ngene     = size(GenesData,1);

%% Pre-evaluate single gene classifiers
% Here I evaluate all possble single-gene classifiers (AtomicNets)
PrimeNets.Function  = Construct_Function([1:Ngene -(1:Ngene)]);
PrimeNets.IsaNOT    = [false(Ngene,1); true(Ngene,1)]; % Just denoting that this is a "NOT" circuit
PrimeNets.GeneID    = [(1:Ngene)'; (1:Ngene)'];

% Remove genes with negligible effect
ExpRanges(:,1) = min(GenesData,[],2);
ExpRanges(:,2) = max(GenesData,[],2);
tmpReponse = log2(Evaluate_Function_C(PrimeNets.Function, ExpRanges));
ResponseRange = range(tmpReponse,2); % this is the most effect this single-gene function can make alone
PrimeNets = Pej_Struct_RowSelect(PrimeNets, ResponseRange>Consts.Learning_mRT_log2); clear NFhits

% the rest
AtNOutput    = log10(Evaluate_Function_C(PrimeNets.Function, GenesData));
Best_Healthy = min(AtNOutput(:,~Annotation),[],2);  % This is the bests Healthy output by this gene. 
Best_Cancer  = max(AtNOutput(:, Annotation),[],2);  % This is the bests Cancer output by this gene. 

AtNet_perf = nan(size(AtNOutput)); % This variable says how does this specific circuit classifies this specic sample as compared to the easiest sample in the other class. So basically if this is negative it's hopeless.
AtNet_perf(:, Annotation)   =  AtNOutput(:, Annotation) - repmat(Best_Healthy,1,sum( Annotation));
AtNet_perf(:,~Annotation)   = -AtNOutput(:,~Annotation) + repmat(Best_Cancer ,1,sum(~Annotation));  
AtNet_perf(AtNet_perf<0) = 0;

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


E_Pos = mean(AtNOutput(:, Annotation),2);
E_Neg = mean(AtNOutput(:,~Annotation),2);

Clasification_pwr                = zeros(size(AtNOutput));
Clasification_pwr(:, Annotation) =  AtNOutput(:, Annotation) - repmat(E_Neg, 1, sum( Annotation));
Clasification_pwr(:,~Annotation) = -AtNOutput(:,~Annotation) + repmat(E_Pos, 1, sum(~Annotation));

PrimeNets.Clasification_ub = max(Clasification_pwr,0);
end




