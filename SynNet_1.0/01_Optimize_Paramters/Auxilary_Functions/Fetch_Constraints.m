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
InBuff = textscan(F, '%s%s', 'CommentStyle', '#', 'delimiter', '\t');

if ~isempty(ConstraintName)
    C = Cast_String(InBuff{2}{strcmp( InBuff{1}, ConstraintName)});
    fclose(F);
    return
else
    for i = 1:length(InBuff{1});
        C.(InBuff{1}{i})= Cast_String(InBuff{2}{i});
    end
    fclose(F);
    return
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

