% This file reads a tab-separated(other delimiters need to be given explicitely) meta-data table with one header line containing the name of the metadata
% Lines starting with '#' are disregarded.

%Pej 2014
%------
function MetaData = Pej_Read_MetaData(File, Delimiter)
if nargin < 2; Delimiter= '\t';end

Fin = fopen(File, 'r');
tmpLine = '#';
while tmpLine(1)=='#';tmpLine = fgetl(Fin);end % Skip header lines
Header = regexp(tmpLine, ['([^' Delimiter ']*)'], 'Tokens');
MetaData.Values = textscan(Fin, repmat( '%s', 1, length(Header)), 'CollectOutput',1, 'delimiter', Delimiter, 'CommentStyle', '#'); fclose(Fin);
MetaData.Values = MetaData.Values{1};
MetaData.Labels = vertcat(Header{:});
end