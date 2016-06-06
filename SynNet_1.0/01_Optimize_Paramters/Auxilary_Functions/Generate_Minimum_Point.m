% This function finds for each function, the point in the Domain which the
% minimum values falls.

function Xmin = Generate_Minimum_Point()
global Domain  FunctionArray
Xmin = nan(max(abs(FunctionArray(:))), size(FunctionArray,1));
for f = 1:size(FunctionArray,1)
    CurrFunction = zeros(size(FunctionArray,2), size(FunctionArray,3));
    CurrFunction(:,:) = FunctionArray(f, :,:); % Note of stupidity: Functioninstance is defined as global in both IdealFun(X) and Evaluate_FunctionCnt_Minimal(X, theta)
    Fun_Genes = unique(CurrFunction(CurrFunction(:)~=0));  % get all the genes, I assume a gene is not used twice
    if length(Fun_Genes)~= length(unique(abs(Fun_Genes)))
        warning('Function has a gene in OR and NOT clause at the same time, this procedure is not designed for it!!')
    end
    for i = 1:length(Fun_Genes)
        if Fun_Genes(i)> 0
            Xmin(Fun_Genes(i),f) = Domain(Fun_Genes(i),1); % Put the lowest value for OR clause
        elseif Fun_Genes(i) < 0
            Xmin(-Fun_Genes(i),f) = Domain(-Fun_Genes(i),2); % Put the highest value for NOT clause
        end
        
    end
end

Xmin = 10.^Xmin;
end