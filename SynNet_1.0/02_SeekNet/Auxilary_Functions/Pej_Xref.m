% This file gets a cross-reference, Identifies the original ID and returns
% the IDs for all the other given types in the reference in output. The
% unfound IDs are returned as the original ID. and the orders and the size
% of the input is preserved.

% Pej 2014 March, pejman.m@gmail.com
%--------------
function Ynonunique = Pej_Xref(Xnonunique,Xref_File, LeaveEmpty)
if nargin<2 || isempty(Xref_File)
    Xref_File = '/Users/pejmanm/Desktop/LocalTMP/PEJ_Resources/Gene-Xref-25Apr2015.txt';
end

if nargin <3
    LeaveEmpty = false; % then it will put the original IDs in the places that they were not found in Xref
end

[X,~,iX] = unique(Xnonunique);
% X = upper(X);
n = length(X);
Xref = Pej_Read_Table(Xref_File, [], false);
IDnames = fieldnames(Xref);
for i = length(IDnames):-1:1
    %     Xref.(IDnames{i}) = upper(Xref.(IDnames{i}));
    Intersct(i) = length(intersect(X(1:min(500,n)), Xref.(IDnames{i})));
end
[nIntersct dbI] = max(Intersct);
disp(['Identified ID Type: ' IDnames{dbI} ' (IDs matched: ' int2str(nIntersct) '/' int2str(min(500,n)) ')']);
[~, ai, bi] = intersect(X, Xref.(IDnames{dbI}));
disp([int2str(length(ai)) ' out of ' int2str(n) ' unique IDs were found in the cross-reference.']);
if length(ai)/n < .9
    beep
    warning('Less that 90% of the IDs were found in the Xref, you might need a better reference!!')
end
for i = length(IDnames):-1:1
    if i~= dbI
        if LeaveEmpty
            Y.(IDnames{i}) = cell(size(X));  % leave empty those that were not found with the original ID
        else
            Y.(IDnames{i}) = X;%cell(size(X));  % fill those that were not found with the original ID
        end
        Y.(IDnames{i})(ai) = Xref.(IDnames{i})(bi);
    end
end

Ynonunique = Pej_Struct_RowSelect(Y, iX);

end