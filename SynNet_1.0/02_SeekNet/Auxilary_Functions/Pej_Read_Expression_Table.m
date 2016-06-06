% This file reads an expression tab-separated(other delimiters need to be given explicitely) table with one header line and first column
% containing gene names. 
% Lines starting with '#' are disregarded.

%Pej 2014
%------
function GeneData = Pej_Read_Expression_Table(File, Delimiter)
if nargin < 2; Delimiter= '\t';end

Fin = fopen(File, 'r');
tmpLine = '#';
while tmpLine(1)=='#';tmpLine = fgetl(Fin);end % Skip header lines
Header = regexp(tmpLine, ['([^' Delimiter ']*)'], 'Tokens');
InBuff = textscan(Fin, ['%s' repmat( '%f', 1, length(Header)-1)], 'CollectOutput',1, 'delimiter', Delimiter, 'CommentStyle', '#','treatAsEmpty',{'NA','na'}); fclose(Fin);
GeneData.GeneNames = InBuff{1};
GeneData.Expressions = double(InBuff{2}); clear InBuff
GeneData.SampleLabels = vertcat(Header{2:end});
end