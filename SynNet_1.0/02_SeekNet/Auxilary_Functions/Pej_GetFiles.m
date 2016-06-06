% This function gets a pattern and runs it within the linux "find" command,
% and returns the path of the files fitting to the pattern.
% Pej 2013

function Files = Pej_GetFiles(Pattern)
if nargin == 0
    Pattern = '*.counts.txt';
end

FI = find(Pattern=='/', 1, 'last');
if ~isempty(FI)
    Folder = Pattern(1:(FI-1));
    Pattern(1:FI)=[];
else
    Folder = '.';
end

[v F] = system(['find ' Folder ' -type f -name "' Pattern '"']); clear v
tmpFiles = textscan(strrep(F, ' ', ':'), '%s'); tmpFiles = tmpFiles{1};

for i = length(tmpFiles):-1:1
    Files{i} = strrep(tmpFiles{i,1}, ':', ' ');
end

end
