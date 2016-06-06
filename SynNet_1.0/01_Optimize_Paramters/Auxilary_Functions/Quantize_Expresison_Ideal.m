% Description:
% This function, Quanties the data into 3 values, Zero, One

% Outputs:
% "BinData" is one for all values that are larger or equal to the geometric mean
% of the Zero and one level threshold.

% Written By Pejman, 14April2014, Basel
% Pejman.m@gmail.com
%----------------------------------------------

function [BinData ] = Quantize_Expresison_Ideal(Data, StepAt)
BinData = (Data>=StepAt);
end