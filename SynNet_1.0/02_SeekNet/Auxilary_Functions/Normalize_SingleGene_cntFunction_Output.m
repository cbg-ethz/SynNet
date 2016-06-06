% This function gets a continuous output matrix "X" transforms it to log-spacce and normalizes is
% between zeor and one. Assuming that the value of a single-Input function
% is falling between the borders of the function value for a zero expressed
% miRNA, or the whole miRNA pool size, for a Literal or a negative literal
% circuit.

% Pej 2014 July
%---------------

function Z = Normalize_SingleGene_cntFunction_Output(X)
global Consts
F = Construct_Function([1 -1],Consts.MaxAnd, Consts.MaxOr );
OutputLimits = log(Evaluate_FunctionCnt(F, [0, Consts.Total_miRNAs_perCell]));
Z = (log(X) - min(OutputLimits(:))) ./ range(OutputLimits(:));
end