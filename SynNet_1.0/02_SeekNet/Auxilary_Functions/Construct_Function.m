% This Script makes a new function from the given single "Gene".

% Inputs:
%   Gene: An integer input. "Gene" can be positive or negative, negative means negate literal.
%   If Gene is an array of length K, then the output will be a 3D array of
%   K functions with (FunctionArray(k,:,:) corresponding to the k'th input.
%
%   MaxAnd: Maximum number of AND inputs in the function family (number of rows);
%
%   MaxOr : Maximum number of OR  inputs in the function family (number of columns);
%
% Written By Pejman, 14Dec2012, Engelberg
% Pejman.m@gmail.com
%----------------------------------------------

function FunctionArray = Construct_Function(Gene)
global Consts

K = length(Gene);
if Consts.MaxOr== 1
        FunctionArray = zeros(K,Consts.MaxAnd,2);
else
    FunctionArray = zeros(K,Consts.MaxAnd,Consts.MaxOr);
end
for k = 1:K
    FunctionArray(k,1,1) = Gene(k);
end

Practicals = Ispractical(FunctionArray);
FunctionArray = FunctionArray(Practicals, :,:);
end