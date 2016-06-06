function [Performance, FPrate, FNrate] = Get_Binary_Performance(FunctionArray, GenesData, TrueAnnoation, FigurePath)

FN  = size(FunctionArray,1);
Ann = repmat(TrueAnnoation,FN, 1);

[BinData, BinDataIsVague] = Quantize_Expresison(GenesData);
[BinResult, BinResultIsVague] = Evaluate_Function_B(FunctionArray ,BinData, BinDataIsVague);
MisClasifications = xor(BinResult, Ann);
MisClasifications = MisClasifications | BinResultIsVague; 

% Classification accuracy
Performance = 100 - sum(MisClasifications,2) ./ size(MisClasifications,2) * 100;
% The false positive rate is FP / (FP + TN).
FPrate = sum(MisClasifications & Ann==false,2) / sum(Ann == false)*100;
% The false positive rate is FN / (FN + TP).
FNrate = sum(MisClasifications & Ann==true ,2) / sum(Ann == true )*100;
if nargin > 3
    % Do nothing for the moment
end
end