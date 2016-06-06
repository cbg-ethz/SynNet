% This function concatenates two structrures of the same format by adding
% one under the other one!
% If A field in a structure has multiple dimensions this always adds over the
% first dimension.

% Pej - Apr 2015
%-----------------------------
function JStruct = Pej_Struct_Cat(Struct1, Struct2)
if isempty(Struct1)|| isempty(Struct2)
    if isempty(Struct1)
        warning('Struct1 is empty!')
        JStruct = Struct2;
    else
        warning('Struct2 is empty!')
        JStruct = Struct1;
    end
    return
end

Fs = fieldnames(Struct1);
Fs2 = fieldnames(Struct1);
if ~isequal(sort(Fs), sort(Fs2))
    error('Incompatible structure!')
end
for i = 1:length(Fs)
    JStruct.(Fs{i})= [Struct1.(Fs{i});Struct2.(Fs{i})];
end
end