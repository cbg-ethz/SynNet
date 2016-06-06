function FunctionArray = FixFunctionFormat(FunctionArray, FormattedHit)
if length(size(FunctionArray))==2; FunctionArray = reshape(FunctionArray, [1,size(FunctionArray)]);end
if nargin < 2; FormattedHit = false(size(FunctionArray,1));    end

tmpHitPool = squeeze(FunctionArray(1,:,:));
tmpFunct   = zeros(size(FunctionArray,1),length(tmpHitPool(:)));
for i = size(FunctionArray,1):-1:1
    %% Reduce each function individually
    clear currCF
    currCF(:,:) = FunctionArray(i,:,:);
    if FormattedHit(i)
        % I assume it's well-structured already!
        tmpHitPool = currCF;
    else
        % Reformat, and reduce the function
        currCF_RF = Reformat_Function(currCF);
        currCF    = Reduce_Function(currCF_RF);
        tmpHitPool(:,:)        = 0;
        tmpHitPool(1:size(currCF,1),1:size(currCF,2))     = currCF;
        FunctionArray(i,:,:)      = tmpHitPool;
    end
    tmpFunct(i,:) = tmpHitPool(:);
end

S2 = size(FunctionArray,2); S3 = size(FunctionArray,3);
for i = size(tmpFunct,1):-1:1
    if FormattedHit(i)
        FunctionArray(i,:,:) =   FunctionArray(i,:,:);
    else
        FunctionArray(i,:,:) = reshape(tmpFunct(i,:),S2,S3);
    end
end
FunctionArray(size(tmpFunct,1)+1:end,:,:)=[];
end


function RFmyFcn = Reformat_Function(myFcn)
%% Sort the structures
% Such that
% 1: all-zero rows will be removed
myFcn(sum(myFcn~=0,2)==0,:)=[];
% 2: non-zero elemnts are in the begining of the row.
% 3: The numbers in each row are sorted in ascending order
% Additionally, in case of  P OR  P: leave only one of them.
myFcn(myFcn==0)= NaN;
for r = 1:size(myFcn,1)
    tmpRw = unique(myFcn(r,:));
    myFcn(r,:) = NaN;
    myFcn(r,1:length(tmpRw)) = tmpRw;
end
myFcn(isnan(myFcn))= 0;
% 4: The non-zero rows are sorted in descending order based on the first element in the row.
% Additionally: In case of  P AND  P: leave only one of them.
RFmyFcn = unique(myFcn,'rows');
end

function myFcn = Reduce_Function(myFcn)
%% Eliminate redundancies in a formatted function
% In case of ~P AND P: delete function by setting all to zero.
% WARNING: This part is only correct under the assumption that NOT gate
% does not accur within and OR gate, i.e. if there's a NOT the rest of the
% columns are assumed to be empty.

if  ~isequal(size(unique(abs(myFcn),'rows'),1), size(myFcn,1))
    myFcn(:) = 0;
    return
end


% In case of P AND (P OR Q): just keep the P
R = size(myFcn,1);
if R<2; return;end;
for r1 = 1:R
    row1 = myFcn(r1,myFcn(r1,:)~=0);
    LR1  = length(row1);
   %  warning('CHECK HERE')
    if LR1 == 1
        myFcn(-myFcn==row1) = 0; % Not P AND  (P OR Q) -> Not P AND Q
    end
end
for r1 = 1:R
    row1 = myFcn(r1,myFcn(r1,:)~=0);
    LR1  = length(row1);
    for r2 = r1+1:R
        row2 = myFcn(r2,myFcn(r2,:)~=0);
        if LR1 + length(row2) > 2
            r1r2Intsct = intersect(row1, row2);
            if isequal(size(row1), size(r1r2Intsct)) || isequal(size(row2), size(r1r2Intsct))
                % Then one is a subset of another, just keep the intersect.
                myFcn(r1,:) = 0;
                myFcn(r1,1:length(r1r2Intsct)) = r1r2Intsct;
                row1 = r1r2Intsct;
                
                myFcn(r2,:) = 0;
            end
        end
    end
end
end



