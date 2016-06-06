% This code is a simple form of the matlab 'intersect' command.
% The main difference is that it assumes the input arrays are pre-sorted in
% Ascending order. This is much fasted the "intersect" in repetitive tasks


% Pej Apri 2015, NYGC
%------------------------

function [IDX] = Pej_Intersect_SortedVectors(X,Y)

i1 = 1;i2=1;
l1 = length(X);l2 = length(Y);

IDX = nan(size(X));

if isnumeric(X)
    while i1<=l1 && i2<=l2
        d=X(i1)-Y(i2);
        
        if d==0
            IDX(i1)=i2;
            i1=i1+1;
            i2=i2+1;
        elseif d<0
            i1=i1+1;
        else % d>0
            i2=i2+1;
        end
    end
else
    while i1<=l1 && i2<=l2
%         [d]=Pej_strcmp(X{i1},Y{i2});
        [d]=strcmpC(X{i1},Y{i2});
        
        if d==0
            IDX(i1)=i2;
            i1=i1+1;
            i2=i2+1;
        elseif d<0
            i1=i1+1;
        else % d>0
            i2=i2+1;
        end
    end
    
end

IDX(isnan(IDX))=[];


