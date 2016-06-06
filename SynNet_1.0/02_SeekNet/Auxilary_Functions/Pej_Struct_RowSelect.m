% This function selects rows, from fields in a structure
% It's assumed that all fields have the same  number oof rows.
% "IndextoSelect" is a binary vector of the same size as the filed, or the
% indices of the rows to be selected.
% Example:
%
% X  = Pej_Read_Table('A_Tab-separated_file.txt');
% X = Pej_Struct_RowSelect(X, [1 5 12]) % To Select 1st, 5th and 12th rows
% from all fileds,

% Pej - Oct 2013
function Struct = Pej_Struct_RowSelect(Struct, IndextoSelect)

Fs = fieldnames(Struct);
for i = 1:length(Fs)
    Sz = size(Struct.(Fs{i}));
    switch length(Sz)
        case 1
            Struct.(Fs{i}) = Struct.(Fs{i})(IndextoSelect);
        case 2
            Struct.(Fs{i}) = Struct.(Fs{i})(IndextoSelect,:);
        case 3
            Struct.(Fs{i}) = Struct.(Fs{i})(IndextoSelect,:,:);
        otherwise
            Szn = [sum(IndextoSelect>0), Sz(2:end)]; % Update the first dimension size, and reshape it to the original dimensionality.
            Struct.(Fs{i}) = reshape(Struct.(Fs{i})(IndextoSelect,:),Szn);
    end
end
end