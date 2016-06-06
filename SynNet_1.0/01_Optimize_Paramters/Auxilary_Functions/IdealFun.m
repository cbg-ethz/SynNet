function Y = IdealFun(X)
%% Target function
global Step_High Step_Low oBP
BinData = Quantize_Expresison_Ideal(X, oBP);
% BinOut = Evaluate_Function(Functioninstance ,BinData);
BinOut = Evaluate_Function_Minimal(BinData);
Y = (Step_High-Step_Low) * BinOut+Step_Low;

end