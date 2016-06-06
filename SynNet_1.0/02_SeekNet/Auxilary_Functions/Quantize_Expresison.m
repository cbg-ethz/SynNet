% Description:
% This function, Quanties the data into 3 values, Zero, One, and NaN(Vague)
% This is coded by two Matrices, one binary matrix, "BinData", of values, and one extra
% "IsVague" binary matrix, that denotes the vague values.

% Outputs:
% "BinData" is one for all values that are larger or equal to the geometric mean
% of the Zero and one level threshold.
% "IsVague" is true, for the points that are between "ZeroLvl" and "OneLvl"

% Written By Pejman, 14Dec2012, Engelberg
% Pejman.m@gmail.com
%----------------------------------------------

function [BinData, IsVague, Margin, MidThr] = Quantize_Expresison(Data)
global Consts

MidThr  = exp(mean([log(Consts.Quantizer_ZeroLvl), log(Consts.Quantizer_OneLvl)])); % Geo-mean of the thrholds
BinData = (Data>=MidThr);
IsVague = (Data>Consts.Quantizer_ZeroLvl) & (Data<Consts.Quantizer_OneLvl);

Margin  = nan(size(Data));
Margin( BinData)  =  log10(Data( BinData)) - log10(MidThr);
Margin(~BinData)  = -log10(Data(~BinData)) + log10(MidThr);
end