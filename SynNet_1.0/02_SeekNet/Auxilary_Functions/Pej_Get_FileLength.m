% This code is equal to "wc -l TextFile" in linux/unix
% Pej Now 2013
% -----
function FLength = Pej_Get_FileLength(TextFile)

fid = fopen(TextFile, 'r');
[~,  FLength] = fscanf(fid,'%*[^\n]%1c');

end