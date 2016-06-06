% This function gets an array of classifier functions (See
% Evaluate_Function.m and Construct_Function.m for the format details and for deffinition of the classifier function),
% and:
% 1: returns a boolean array, "Practicals",  indicating the paractical functions
% under the secondary constraints indicated in the global
% variable "Consts", these secondary constraints include MaxOrCount, MaxNotCount
% and MaxCircuitSize, and the fact that NOT gate and OR gates cannot be
% combined.

% The function works in two "Mode"s: 'Check'(Default), and 'Prune'. The
% Check mode simply Checks the invalid hits, in this mode, the second
% output is empty.
% in Prune mode, function removes invalid parts from the function till it's
% valid. in this mode the second output includes the Pruned funtion array.

% Written by pejman
% 24th of Oct 2013 Basel/July15th2014
% ------------------------------

function [Practicals, FunctionArray] = Ispractical(FunctionArray, Mode)
global Consts
if nargin<2; Mode   = 'Discard'; end
if length(size(FunctionArray))==2; FunctionArray = reshape(FunctionArray, [1,size(FunctionArray)]);end

N     = size(FunctionArray,1);
sFn  = sum(sign(FunctionArray),3);
sFn2 = sum(FunctionArray~=0   ,3);


%% Not "NOT" and "OR" gate mixed, and There's not more than one NOT gate in a row,
RowValid  = (abs(sFn)==sFn2) & (sFn>=-1);
structValid  = all(RowValid,2);
%% MaxOrCount
orValid  = sum(sFn>0,2)<= Consts.MaxOrCount;

%% MaxNotCount
notValid = sum(sFn<0,2)<= Consts.MaxNotCount;

%% MaxCircuitSize
sizeValid = sum(abs(sFn),2)<= Consts.MaxCircuitSize;

Practicals = structValid & orValid & notValid & sizeValid;

%% Prune if Necessary
if strcmpi(Mode, 'Discard') || all(Practicals)
    % Do Nothing
elseif strcmpi(Mode, 'Prune')
    InvalidI(1,:) = find(~Practicals);
    for f = InvalidI
        % Remove one of the problems
        if ~structValid(f)
            FunctionArray(f,~RowValid(f,:),:) = 0; % remove the invalid branch(es) from the AND
        elseif ~sizeValid(f)
            tmpF = find(FunctionArray(f,:));
            FunctionArray(f,randperm(length(tmpF),length(tmpF)-Consts.MaxCircuitSize)) = 0; % remove randomly the extra genes
        else
            if ~orValid(f)
                tmpF = find(sFn(f,:)>0);
                FunctionArray(f,tmpF(randperm(length(tmpF),length(tmpF)-Consts.MaxOrCount)),:) = 0; % remove randomly the extra ORs
            end
            if ~notValid(f)
                tmpF = find(sFn(f,:)<0);
                FunctionArray(f,tmpF(randperm(length(tmpF),length(tmpF)-Consts.MaxNotCount)),:) = 0; % remove randomly the extra NOTs
            end
        end
    end
    try
    [~, FunctionArray(InvalidI,:,:)] = Ispractical(FunctionArray(InvalidI,:,:),'Prune');
    catch
        fprintf('Ispractical crashed:\n%s\n', SPrint_Function(FunctionArray(InvalidI,:,:)));
        FunctionArray(InvalidI,:,:) = 0;
       
    end
else
    error('Invalid Mode of function! Valid Values: [Discard/Prune]')
end
end