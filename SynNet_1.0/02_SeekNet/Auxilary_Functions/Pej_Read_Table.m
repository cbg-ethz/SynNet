% This reads a .CSV or a tab-separated table.
% The first line is expected to be the header
% Written by PEj oct 2013
%-------------

function Data = Pej_Read_Table(Path2File, DLM, ReportSwitch, FormatString)

Fin = fopen(Path2File, 'r');
% Identify Names
RawHeader = fgetl(Fin);
if nargin==1 || isempty(DLM)
    DLM = '[,\t]'; % list of potential delimiters
    dlm    = regexp(RawHeader, DLM, 'match', 'once');
else
    dlm = DLM;
end

if nargin<3 || isempty(ReportSwitch)
    ReportSwitch= true;
end

Header = regexp(RawHeader, dlm, 'split');
Header = Remove_nonAlphanumerics(Header);
noHeaderLines = 1; % Number of header lines
%% Identify File format
if nargin<4 || isempty(FormatString)
    %% Identify the format
    Firstrow = regexp(fgetl(Fin), dlm, 'split');
    % fclose(Fin);
    if Firstrow{1}(1)=='#'
        % This is a header row with values to guide correct format specification.
        % For example if the fisrt row of a numerical column by chance is "NaN",
        % auto formating would recognize it as "Text", this extra line that startes
        % with "#' holds values that ensure correct format recognition, without
        % being included in the final data.
        Firstrow{1}(1)=[];
        noHeaderLines = 2;
    end
    FormatS = '';
    Report  = '';
    if length(Header)<length(Firstrow)
        warning('The header line does not match the data')
        for c = length(Header)+1:length(Firstrow)
            Header{c} = ['UnLabeled_C' int2str(c)];
        end
    end
    for i = 1:length(Firstrow)
        while strcmpi('NA', Firstrow{i}) || strcmpi('NaN', Firstrow{i})
                        % it might be a missing numerical or missing text
            Firstrow = regexp(fgetl(Fin), dlm, 'split');
        end
        %     if ~isempty(regexp(Firstrow{i}, '^[+-]?\d+$'))
        %         % It's an integer
        %         FormatS = [FormatS '%f' ];
        %         Report  = [Report  Header{i} '(int)\t'];
        %     elseif ~isempty(regexp(Firstrow{i}, '^[+-]?\d+(\.)?\d*[eE]?[+-]?\d+$'))
        %         %read it as Float
        %         FormatS = [FormatS '%f' ];
        %         Report  = [Report  Header{i} '(float)\t'];
        %     else
        %         %read it as string
        %         FormatS = [FormatS '%s' ];
        %         Report  = [Report  Header{i} '(str)\t'];
        %     end
        if ~isnan(str2double(Firstrow{i}))
            %read it as Float
            FormatS = [FormatS '%f' ];
            Report  = [Report  Header{i} '(numeric)\t'];
        else
            %read it as string
            FormatS = [FormatS '%s' ];
            Report  = [Report  Header{i} '(str)\t'];
        end
    end
    fclose(Fin);
    if ReportSwitch
        fprintf(['Input file: %s \nRecognized format:\n' Report '\n'], Path2File);
    end
else
    %% Use the provided format
    FormatS=FormatString;
end
%% Read file
Fin = fopen(Path2File, 'r');
DES_Res = textscan(Fin, FormatS, 'HeaderLines',noHeaderLines,'Delimiter', dlm, 'treatAsEmpty',{'NA','na'});%, 'bufsize', 1E+5);
for i = 1:length(DES_Res)
    Data.(Header{i}) = DES_Res{i}; DES_Res{i} = [];
end

%% Make sure the whole file was read thorugh
if ~feof(Fin)
    beep
    warning('Import stopped before reaching the end of the file. Here''s the first unread lines:\n%s\n%s\n', fgetl(Fin), fgetl(Fin));
end
fclose(Fin);
end


function S = Remove_nonAlphanumerics(S)
% this removed all but alphanumerics and substitutes them with '_'. Extra
% underscores at the begining and the end are removed.
for i = 1:length(S)
    tmpS = S{i};
    Valid = ...
        (tmpS >='0' & tmpS <= '9')  | ...
        (tmpS >='a' & tmpS <= 'z')  | ...
        (tmpS >='A' & tmpS <= 'Z')  | ...
        (tmpS == '_');
    tmpS(~Valid)= '_';
    while ~isempty(tmpS)&& tmpS(1)  =='_'; tmpS(1)  =[];end
    while ~isempty(tmpS)&&tmpS(end)=='_'; tmpS(end)=[];end
    if isempty(tmpS); tmpS=['UnLabeled_C' int2str(i)];end
    S{i} = tmpS;
end

for i = 1:length(S)
    if S{i}(1)>='0' && S{i}(1)<='9'
        % There's a digit at the begining
        S{i} = [ 'F_' S{i}];
    end
end
end
