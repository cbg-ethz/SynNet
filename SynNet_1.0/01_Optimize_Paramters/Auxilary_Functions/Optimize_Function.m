function Problem = Optimize_Function(Problem)
set(0, 'defaultTextInterpreter', 'tex');

global Step_High Step_Low Pa Pb Pc PDname Domain optTol FunctionArray FunctionWeight 
global Pscalefactor PosScale InitialParams MaxParamRange InitCost Class_Mass_Balance
Domain = Problem.Log10Expressions.Domain;
Step_High = Problem.IdealFunction.Step_High;
Step_Low  = Problem.IdealFunction.Step_Low;
Pa = Problem.Log10Expressions.Pa;% First  parameter of the Background log10-expression pdf
Pb = Problem.Log10Expressions.Pb;% Second parameter of the Background log10-expression pdf
Pc = Problem.Log10Expressions.Pc;% Third  parameter of the Background log10-expression pdf
PDname = Problem.Log10Expressions.PDF;% Distribution name of log10-expression of the data
FunctionArray  = Problem.IdealFunction.FunctionArray;
FunctionWeight = Problem.IdealFunction.FunctionWeight/sum(Problem.IdealFunction.FunctionWeight);
optTol = log10(Problem.Optimization.Relative_Percision);
stepAt = Problem.IdealFunction.stepAt;
InitialParams = Problem.Continuous_Function.InitialParams;
Class_Mass_Balance = Problem.Optimization.Class_Mass_Balance;
for d = 1:length(Pa)
    % Claculate the scale factor to normalize the truncated probability
    Pscalefactor(d) = 1 / (cdf(PDname{d},Domain(d,2),Pa(d), Pb(d), Pc(d)) -cdf(PDname{d},Domain(d,1),Pa(d), Pb(d), Pc(d)));
end
% Initialize the relative mass of positive/negative samples for each function
PosScale = nan(size(FunctionArray,1),1);
InitCost = nan;
try
    % New MAtlab versions
    OPToptions = optimoptions(@fmincon, 'MaxFunEvals', 3000, 'MaxIter', 1000, 'GradObj', 'off', 'UseParallel', never ,  ...
        'TolX', 10^(optTol) , 'TolFun', 10^(optTol));%, 'Diagnostics', 'on', 'Display','iter-detailed', 'PlotFcns', @optimplotx);
catch
    % Older Matlab versions
    OPToptions =     optimset('MaxFunEvals', 3000, 'MaxIter', 1000, 'GradObj', 'off', 'UseParallel', 'never' , ...
        'TolX', 10^(optTol) , 'TolFun', 10^(optTol), 'Diagnostics', 'off', 'Display','iter-detailed', 'PlotFcns', @optimplotx);
end

Problem.objective = @(x)Cost_function(stepAt, x);
Problem.x0 = (InitialParams);
Problem.solver = 'fmincon';
Problem.options= OPToptions;
Problem.lb = Problem.Continuous_Function.Params_LowerBound;
Problem.ub = Problem.Continuous_Function.Params_UpperBound;
tic
Mabbr = max([InitialParams - Problem.lb, Problem.ub - InitialParams],[],2);% Maximum abbaration from the Initial values
MaxParamRange = sum(Mabbr.^2); 

[Problem.OPTtheta,~,Problem.exitflag] = fmincon(Problem);
Problem.OPTtheta
toc
%% Plotting
Visualize_Function(FunctionArray, Problem)
end

