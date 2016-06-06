% This function removes identical functions from the pool. This is a
% special case for the more general function  "Refine_HitPool"
% NOTE: This assumes all the pool is already formatted and also is sorted
% by performace (if the 3rd input is missing), so it only checks if consecutive functions are equal as
% they are. If 3rd input is there, the function sorts by "Sorting_Measure",
% and then looks for consecutive repeats with equal sorting measures.

function Uniques = Remove_ConsecutiveRepeats(HitPool, Sorting_Measure)
N = size(HitPool,1);
Uniques = true(N,1);

if nargin >1
    if length(Sorting_Measure)~=N; error('Sorting measure does not fit to hitpool!'); end
    % Use the sorting measure, to sort
    [S, Ix] = sort(Sorting_Measure, 'descend');
    Sd0 = (S(1:N-1)-S(2:N))==0;
    
    for i = 1: N-1
        if Sd0(i)
            % Then "i"th and "i+1"th entry have the same value for the sorting measyure, check
            % if they are equal
            CurrFun = HitPool(Ix(i),:,:);
            j=i+1;
            while j<= N && Sd0(j-1)
                t = Ix(j);
                % Look in all function with equal "Sorting measure". This
                % is important because two different function might hapen
                % to give eual results and sort after each other.
                if Uniques(t)
                    if isequaln(CurrFun, HitPool(t,:,:))
                        Uniques(t)=false;
                    end
                end
                j=j+1;
            end
        end
    end
else
    
    for i = 1: N-1
        if isequaln(HitPool(i,:,:), HitPool(i+1,:,:))
            Uniques(i+1)=false;
        end
    end
end
end