% this function reads in the constraints from a text file.
% Usage:

% 1st Form:
% C = Fetch_Constraints(Constraints_File); Fetches all the constraints in
% the files specified in the TAB-separated file "Constraints_File" and puts
% them as separated field in C. The file should be formatted like this:

% Const1    Const1Value
% Const2    Const2Value
% # Comment like, blahblah
% Const3    Const3Val

% the file does not need to be ordered.

% 2nd Form:
% ConstraintValue = Fetch_Constraints(Constraints_File, ConstraintName); Fetches the
% constraint specified by "ConstraintName"

% Written by Pejman
% Oct 22nd, 2013
function C = Fetch_Constraints(Constraints_File, ConstraintName)
if nargin == 0; Constraints_File = 'Constraints.txt';end
if nargin  < 2; ConstraintName = [];end

F = fopen(Constraints_File, 'r');
InBuff = textscan(F, '%s%s', 'CommentStyle', '#', 'delimiter', '\t', 'MultipleDelimsAsOne', 0);

if ~isempty(ConstraintName)
    % The mode that function is called to get a specific contrant
    C = Cast_String(InBuff{2}{strcmp( InBuff{1}, ConstraintName)});
    fclose(F);
else
    % The mode that function is called to get all of the constraints in the
    % input file
    for i = 1:length(InBuff{1});
        C.(InBuff{1}{i})= Cast_String(InBuff{2}{i});
    end
    fclose(F);
    %% Put default values
    C.AnalysisMode = upper(C.AnalysisMode);
    if ~strcmp(C.AnalysisMode,'C') && ~strcmp(C.AnalysisMode,'B')
        error('Invalid Constraint value: "AnalysisMode" value can only be "C"(Continuous) or "B"(Binary).'); 
    end
    
    if ~isfield(C, 'Total_miRNAs_perCell') || isempty(C.Total_miRNAs_perCell)
        C.Total_miRNAs_perCell = 25000;
    end
    
    if ~isfield(C, 'Learning_Avrage_Margin_Weight') 
        C.Learning_Avrage_Margin_Weight = .5;
    end
    
    if ~isfield(C, 'Learning_Worst_Margin_Weight') 
        C.Learning_Worst_Margin_Weight = .5;
    end
    
    if ~isfield(C, 'Learning_Convergence_Thr_2')
        C.Learning_Convergence_Thr_2 = C.Learning_Convergence_Thr_1;
    end
    
    if ~isfield(C, 'Learning_mRT_log2')
        C.Learning_mRT_log2 = .1;
    end
    
    if ~isfield(C, 'Learning_Pruning')
        C.Learning_Pruning = 'off';
    end
    C.Learning_Pruning = lower(C.Learning_Pruning);
    
[~,tmpOutFolder ,~] = fileparts(Constraints_File);
C.OutFolder =  ['F03_Results/' tmpOutFolder '_Results'];
       
    if ~isfield(C, 'DataPath') || isempty(C.DataPath)
        % Its a simulation
        %C.OutFolder =  ['F03_Results/Simulation_' C.AnalysisMode]; 
        C.IsSimulation = true;
        if ~isfield(C, 'Simulation_DataSynthesis_ratio')
            C.Simulation_DataSynthesis_ratio=1;
        end
        
    else
        % it's real data analysis
        C.IsSimulation = false;
        C.Simulation_OrderOnly = 0;
        %[~,tmpOutFolder ,~] = fileparts(Constraints_File);
       % C.OutFolder =  ['F03_Results/' tmpOutFolder '_' C.AnalysisMode];
        C.JointReportName = [C.OutFolder '/' 'Summary_Report_' tmpOutFolder '.txt'];
    end
    
    C.Learning_mRT_log2 = max(eps, C.Learning_mRT_log2);
    
    C.Quantizer_OneLvl  = C.Quantizer_Threshold;
    C.Quantizer_ZeroLvl = C.Quantizer_Threshold-eps;
    
    C.Plotting = ~strcmpi('off', C.Plotting);
    
    Total_MW = (C.Learning_Worst_Margin_Weight + C.Learning_Avrage_Margin_Weight);
    C.Learning_Worst_Margin_Weight  = C.Learning_Worst_Margin_Weight  / Total_MW;
    C.Learning_Avrage_Margin_Weight = C.Learning_Avrage_Margin_Weight / Total_MW;
end


end


function T = Cast_String(S)
% This function gets a cell and converts it into string, or to a number if
% it's numerical.
if nargin ==0; T=[]; return; end
if isempty(S); T=S ; return; end
[T  Status] = str2num(S);
if Status == true
    % It's ok, it was a number
    return
else
    % leave it as it is
    T = S;
end
end

