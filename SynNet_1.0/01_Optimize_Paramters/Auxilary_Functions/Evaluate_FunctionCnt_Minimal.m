% This file evaluates the continous for of a boolean function over values in "X",
% This function is supposed to be faster than the full version.

% Params = [
%     Consts.Continuous_Circuit_F1C
%     Consts.Continuous_Circuit_F2C
%     Consts.Continuous_Circuit_Tmax
%     Consts.Continuous_Circuit_FF4max
%     Consts.Continuous_Circuit_OUTmax];

% Written By Pejman, April2014, Basel
% Pejman.m@gmail.com
%----------------------------------------------

function Result= Evaluate_FunctionCnt_Minimal(X, Params)
global Functioninstance

MaxAnd = size(Functioninstance,1); % Number of AND inputs in the function family
N      = size(X,2);             % Number of Samples in the data.

RowValues  = zeros(MaxAnd, N);
Fneqs = Functioninstance(:)<0;% The NOT part
%% handle the OR part
for i = 1:MaxAnd
    Filt = Functioninstance(i,:)>0;
    if any(Filt)
        RowValues(i,:)= F1(X(Functioninstance(i, Filt),:),Params(1));
    end
end
mirFF4 = Params(4)* F2(sum(RowValues,1), Params(2), Params(3));
%% handle the NOT part
Result = Params(5) * F1([mirFF4;X(-Functioninstance(Fneqs),:)],Params(1));
end


function Y = F1(X,Continuous_Circuit_F1C)
S = sum(X,1);
Y = 1 - S./(S+Continuous_Circuit_F1C);
end

function Y = F2(X, Continuous_Circuit_F2C, Continuous_Circuit_Tmax)
Y = X./(X+Continuous_Circuit_F2C/Continuous_Circuit_Tmax);
end