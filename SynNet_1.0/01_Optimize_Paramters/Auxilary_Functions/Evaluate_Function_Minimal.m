% This file evaluates the boolean function over values in "X", and if
% provided "IsVagueX".

% Inputs:
%   FunctionArray: A single boolean function or an array of K functions with each FunctionArray(k,:,:) corresponding to
%   one signle function. The function is encoded like:
%   F = [5 6 0 0 ; -3 0 0 0 ; 0 0 0 0 ; 0 0 0 0]
%   with each row corresponding inputs of one OR, and all the OR outputs then connected with an
%   AND together, and zeros denoting the unused inputs. Hence, the above example
%   would be interpreted as:
%   F = (Gene5 OR Gene6) AND (NOT Gene3)
%   Additionally we see that F is from a family of functions with maximum 4 inputs in ORs(number of columns)
%   and 4 inputs to the ANDs(number of rows).
%
%   X: This is the gene data matrix. Each row corresponds to one gene, and
%   each column is one sample. Matrix is assumed to be binary.
%
%   IsVagueX(Optional): This is a binary matrix, the same size as "X",
%   denoting the vague values in the "X". For example if "X" has been quantized
%   by X(X<100)=0 and X(X>1000)=1, then all values between 100 and 1000 are vague.
%   here we try to to evaluate the function by disregarding these vague
%   values, if possible. If not, the final output would be also Vague.

% Otputs:
%   Result: A K times N binary matrix reporting the F(X) with each row corresponding to
%   FunctionArray(k,:,:) and each column corrsponding to a sample(i.e. X(:,n)).
%
%   IsVague: Binary matrix of the same size as "Result". If the result of the function is vague.
%   In this case the Result value should be disregarded.
%
%   RowValues: A MaxAnd times N binary matrix, with MaxAnd=number of AND inputs in the function family.
%   Presents the values for each individual OR, ONLY for the first function
%   (i.e. FunctionArray(1,:,:))
%
%   IsVagueRow: Like "IsVague", but for the "RowValues".

% Written By Pejman, 14Dec2012, Engelberg
% Pejman.m@gmail.com
%----------------------------------------------

function Result = Evaluate_Function_Minimal(X)
global Functioninstance


MaxAnd = size(Functioninstance,1); % Number of AND inputs in the function family
N      = size(X,2);             % Number of Samples in the data.

%% Old evaluator
% RowValues  = false(MaxAnd, N);
% for i = 1:MaxAnd
%     Filt = find(Functioninstance(i,:));
%     if isempty(Filt)
%         RowValues(i,:)= true;
%     else
%         for j = Filt
%             xF = Functioninstance(i,j);
%             if xF % If it is not zero, i.e. if it's connected to something
%                 tmpFeat = false(1,N);
%                 if xF<0
%                     % Negated literal
%                     tmpFeat = ~X(-xF,:);
%                 else
%                     tmpFeat =  X( xF,:);
%                 end
%                 RowValues(i,:)  = RowValues(i,:)|tmpFeat;
%             end
%         end
%         % OR lines that are zero when Vague inputs are excluded, will remain vague
%         % So I put one in them, and see if the downstream AND is already
%         % flase without this, otherwise the final output of the AND will
%         % also be vague.
%     end
%     
% end
% Result_old = all(RowValues,1);
% 

%% New evaluator
% should be faster but the same
RowValues  = true(MaxAnd, N);
Fpos = max(0,Functioninstance); % The OR  part
Fneqs = Functioninstance(:)<0;  % The NOT part

%% handle the OR part
BranchCnt = 0;
for i = 1:MaxAnd
    Filt = Fpos(i,:)>0;
    if any(Filt)
        % It's not an empty Branch
        BranchCnt = BranchCnt+1;
        RowValues(BranchCnt,:)= any(X(Functioninstance(i, Filt),:),1);
    end
end
%% handle the NOT part
BranchCnt = BranchCnt+1;
RowValues(BranchCnt, :) = all(~X(-Functioninstance(Fneqs),:),1);

Result = all(RowValues(1:BranchCnt,:),1);
% 
% if ~isequal(Result_old, Result)
%     warning('THe new and old code don''t match!!')
% end
% disp('This code is redundant!!')
        
        
end