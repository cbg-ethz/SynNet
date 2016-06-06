% This just reads Bed files.
% Reference: http://genome.ucsc.edu/FAQ/FAQformat.html#format1


% Pejman 2014
%--------------

function Data = Pej_Read_Bed(Path2File)

DLM = '[\t]'; % list of potential delimiters
fprintf(['Input file: %s'], Path2File);

FormatS = '%s%d%d%s%d%c%[^\n]';
Header = {'chr', 'start', 'end', 'name', 'score', 'dir', 'the_rest'};
%% Read file
Fin = fopen(Path2File, 'r');
DES_Res = textscan(Fin, FormatS,'Delimiter', DLM, 'commentStyle','#');%, 'bufsize', 1E+5);
for i = 1:length(DES_Res)
    Data.(Header{i}) = DES_Res{i}; DES_Res{i} = [];
end
fclose(Fin);
fprintf('\tdone!\n')
end



