% This Code optimized the biochemical paramters over the circuits defined
% in "Problem.IdealFunction.FunctionArray" around the binarization threshod
% indicated by "Problem.IdealFunction.stepAt'

% by Pejman Mohammadi
% pejman.m@gmail.com
%-----------

function Launch_Global_Optimization()
addpath('Auxilary_Functions')

%% Initial Circuit Constants
Consts.Continuous_Circuit_F1C=	40;
Consts.Continuous_Circuit_F2C=	2000;
Consts.Continuous_Circuit_Tmax=	50000;
Consts.Continuous_Circuit_FF4max=	300;
Consts.Continuous_Circuit_OUTmax=	50000;

%% Initial Settings and bounds for the Continuous Circuit in Log10
Problem.Continuous_Function.InitialParams = log10([
    Consts.Continuous_Circuit_F1C
    Consts.Continuous_Circuit_F2C
    Consts.Continuous_Circuit_Tmax
    Consts.Continuous_Circuit_FF4max
    Consts.Continuous_Circuit_OUTmax
    ]);

Problem.Continuous_Function.Params_LowerBound = [log10(20) ; Problem.Continuous_Function.InitialParams(2:end) - 1]; 
Problem.Continuous_Function.Params_UpperBound = Problem.Continuous_Function.InitialParams + 1;

%% Description of the Ideal Circuit
Problem.IdealFunction.stepAt    = 25000 * 10^-2; % binarization threshold
Problem.IdealFunction.Step_High = 50000; % Output High value ideal (This value is never used in optimization just used ONLY for visualization)
Problem.IdealFunction.Step_Low  = 500;   % Output Low value ideal (This value is never used in optimization just used ONLY for visualization)

Problem.IdealFunction.FunctionArray(1,:,:) =[ % Circuit: (Gene1)
    1    0    0      ;
    0    0    0      ;
    0    0    0      ;
    0    0    0     ];
% ----- How to Define the circuits  -----
%   FunctionArray: A single boolean function or an array of K functions with each FunctionArray(k,:,:) corresponding to
%   one signle function. The function is encoded like:
%   F = [5 6 0 0 ; -3 0 0 0 ; 0 0 0 0 ; 0 0 0 0]
%   with each row corresponding inputs of one OR, and all the OR outputs then connected with an
%   AND together, and zeros denoting the unused inputs. Hence, the above example
%   would be interpreted as:
%   F = (Gene5 OR Gene6) AND (NOT Gene3)
%   Additionally we see that F is from a family of functions with maximum 4 inputs in ORs(number of columns)
%   and 4 inputs to the ANDs(number of rows).

Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: not(Gene1)
    -1    0    0      ;
    0    0    0      ;
    0    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: (Gene1) OR (Gene2)
    1    2    0      ;
    0    0    0      ;
    0    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: not(Gene1) AND (Gene2)
    -1    0    0      ;
    2    0    0      ;
    0    0    0      ;
    0    0    0     ];

Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: (Gene1) AND (Gene2)
    1    0    0      ;
    2    0    0      ;
    0    0    0      ;
    0    0    0     ];

Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: not(Gene1) AND not(Gene2)
    -1    0    0      ;
    -2    0    0      ;
    0    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: (Gene1) OR (Gene2) OR (Gene3)
    1    2    3      ;
    0    0    0      ;
    0    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: [(Gene1) OR (Gene2)] AND (Gene3)
    1    2    0      ;
    3    0    0      ;
    0    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: (Gene1) AND (Gene2) AND (Gene3)
    1    0    0      ;
    2    0    0      ;
    3    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: not(Gene1) AND [(Gene2) OR (Gene3)]
    -1   0    0      ;
    2    3    0      ;
    0    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: not(Gene1) AND (Gene2) AND (Gene3)
    -1   0    0      ;
    2    0    0      ;
    3    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: not(Gene1) AND not(Gene2) AND (Gene3)
    -1   0    0      ;
    -2   0    0      ;
    3    0    0      ;
    0    0    0     ];
Problem.IdealFunction.FunctionArray(end+1,:,:) =[ % Circuit: not(Gene1) AND not(Gene2) AND not(Gene3)
    -1   0    0      ;
    -2   0    0      ;
    -3   0    0      ;
    0    0    0     ];

Problem.IdealFunction.FunctionWeight = ones(1, size(Problem.IdealFunction.FunctionArray,1));
%% Background miRNA Expression Distribution Assumpptions
Ngenes = max(abs(Problem.IdealFunction.FunctionArray(:))); 
Problem.Log10Expressions.Domain = repmat(log10([1 25000]),  Ngenes,1);
Problem.Log10Expressions.PDF= repmat({'normal'}, Ngenes,1);
Problem.Log10Expressions.Pa = ones(Ngenes,1) .* log10(Problem.IdealFunction.stepAt);% First parameter of the Background log10-expression pdf
Problem.Log10Expressions.Pb = ones(Ngenes,1)* 1;% Second parameter of the Background log10-expression pdf
Problem.Log10Expressions.Pc =  nan(Ngenes,1);% (if applicable) Third parameter of the Background log10-expression pdf

%% Other Optimization Constant
Problem.Optimization.Relative_Percision = 1E-6; %
Problem.Optimization.Class_Mass_Balance = true;  % When this is one the cost function is normalized to weight positive and negative samples equally.
warning off
Problem = Optimize_Function(Problem);
save(['Report_' num2str(round(Problem.IdealFunction.stepAt))] , 'Problem');
warning on
clearvars -global
end