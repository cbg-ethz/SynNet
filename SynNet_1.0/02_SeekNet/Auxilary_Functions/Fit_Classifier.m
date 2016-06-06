function Problem = Fit_Classifier(Problem)
global Consts

if Consts.AnalysisMode == 'B'
    Problem = Fit_Boolean_Classifier_Margin(Problem);
else
    Problem = Fit_Continuous_Classifier(Problem);
end
end