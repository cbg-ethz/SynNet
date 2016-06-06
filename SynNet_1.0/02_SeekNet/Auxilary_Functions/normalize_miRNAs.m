% This file normalizes the miRNA data to sum up to roughly 25000 copies per
% cell.
% WARNING: The code is not general, it's tailored.
% Pejman Dec. 2012, Pejman.m@gmail.com
%----------------------------------
function miRNA_Data = normalize_miRNAs(miRNA_Data, SampleFilt, Consts)

if nargin <2
    SampleFilt = 1:size(miRNA_Data.Expressions,2);
    Consts.Total_miRNAs_perCell = 25000;
    Consts.Expression_Pseudocount = 1;
end
miRNAperCell = Consts.Total_miRNAs_perCell; 

if isempty(SampleFilt)
    SampleFilt = true(1,size(miRNA_Data.Expressions,2));
end
Dnormed = Normalize_Columns(miRNA_Data.Expressions(:,SampleFilt));

NRef = median(sum(Dnormed));
Scale = NRef / miRNAperCell;

miRNA_Data.Dnormed = nan(size(miRNA_Data.Expressions));
miRNA_Data.Dnormed(:,SampleFilt) = Dnormed / Scale;
miRNA_Data.Pseudocount = Consts.Expression_Pseudocount / Scale;
end

function X = Normalize_Columns(X,Offset )
F = sum(X>0,2)>(size(X,2)*.75); % just include those that are present in more than half of the samples
if sum(F) < 50
    warning(sprintf('Low quality data: Only %d genes were included in the library normalization.\n This can happen when some samples have extremely low coverage, try excluding low coverage samples from the data.', sum(F)));
end

tmpXL = log(X(F,:));
XL = log(X);    

if nargin<2
    ref     = median(tmpXL,2);
    Xd      = tmpXL - repmat(ref,1, size(tmpXL,2));
    Filt    = ~any(isnan(Xd),2);
    Offset  = median(Xd(Filt,:));
end
if any(exp(-Offset)>10)
    warning('Low quality data: Some expressions are normalized over 10 folds.');
end

if any(~isfinite(Offset))
    beep
    warning(sprintf('NORMALIZATION FAILED!\nOne or more samples failed in normalization, this happens when more than half of the genes used in normalization are not expressed in a given sample.\n try excluding low coverage samples from the data.'));
    fprintf('Switching to normalization by library size. This can lower the quality of normalization.\n')
    SL = sum(X,1);
    Offset = log(SL) - log(median(SL)); 
end

X = exp(XL - repmat(Offset, size(XL,1),1));
end