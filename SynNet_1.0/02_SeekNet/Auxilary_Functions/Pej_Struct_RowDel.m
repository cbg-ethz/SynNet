% This function removes rows(first dimension), from fields in a structure
% It's assumed that all fields have the same  number oof rows.
% "IndextoRemove" is a binary vector of the same size as the filed, or the
% indices of the rows to be removed.
% Example:
%
% X  = Pej_Read_Table('A_Tab-separated_file.txt');
% X = Pej_Struct_RowDel(X, [1 5 12]) % To remove 1st, 5th and 12th rows
% from all fileds,

% Pej - Oct 2013
function Struct = Pej_Struct_RowDel(Struct, IndextoRemove)
Fs = fieldnames(Struct);
for i = 1:length(Fs)
    Sz = size(Struct.(Fs{i}));
    switch length(Sz)
        case 1
            Struct.(Fs{i})(IndextoRemove) = [];
        case 2
            Struct.(Fs{i})(IndextoRemove,:) = [];
        case 3
            Struct.(Fs{i})(IndextoRemove,:,:) = [];
        otherwise
            Szn = [Sz(1) - sum(IndextoRemove>0), Sz(2:end)]; % Update the first dimension size, and reshape it to the original dimensionality.
            Struct.(Fs{i})(IndextoRemove,:) = [];
            Struct.(Fs{i}) = reshape(Struct.(Fs{i}), Szn);
    end
end
end