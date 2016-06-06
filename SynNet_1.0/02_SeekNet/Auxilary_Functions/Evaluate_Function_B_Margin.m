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

function [Result, Margin] = Evaluate_Function_B_Margin(FunctionArray ,X, Xmargins)
global Consts

if ~islogical(X); warning('Input to the binary function evaluuator is not binary!');end
if nargin <3; Xmargins = zeros(size(X)); end; % If there's no info on the vagueness of dta, assume nothing's vague.
if length(size(FunctionArray))==2
    %% It's a Single function, format it into an 3D array of length one.
    tmpFnc = FunctionArray;
    clear FunctionArray
    FunctionArray(1,:,:) = tmpFnc;
    clear tmpFnc
end
MaxAnd = size(FunctionArray,2); % Number of AND inputs in the function family
N      = size(X,2);             % Number of Samples in the data.

% The two following parts are basically identical, I just duplicate them
% not waste time in cases that we do not care about Margins.
if nargin == 2
    %% Skip Margin Calculations
    for K = size(FunctionArray,1):-1:1
        RowValues  = true(Consts.MaxOrCount+1, N);
        Function(:,:) = (FunctionArray(K,:,:));
        Fpos = max(0,Function); % The OR  part
        Fneqs = Function(:)<0;  % The NOT part
        
        %% handle the NOT part
        RowValues(1, :) = all(~X(-Function(Fneqs),:),1);
        %% handle the OR part
        BranchCnt = 1;
        for i = 1:MaxAnd
            Filt = Fpos(i,:)>0;
            if any(Filt)
                % It's not an empty Branch
                BranchCnt = BranchCnt+1;
                RowValues(BranchCnt,:)= any(X(Function(i, Filt),:),1);
            end
        end
        Result(K,:) = all(RowValues,1);
        Margin = [];
    end
elseif nargin == 3
    %% Add Margin Calculations
    for K = size(FunctionArray,1):-1:1
        RowValues  = true(MaxAnd, N);
        RowMargins =  inf(MaxAnd, N);
        
        Function(:,:) = (FunctionArray(K,:,:));
        Fpos =  max(0,Function); % The OR  part
        Fneg = -min(0,Function); % The NOT part
        if isequal(Fpos, Fneg)
            % It's an empty function
            Result(K,:) =  true(1, N);
            Margin(K,:) = zeros(1, N); 
            continue
        end
        
        for i = 1:MaxAnd
            Filt = Fpos(i,:)>0;
            if any(Filt)
                %% handle the OR part
                % It's not an empty Branch
                RowValues(i,:)    =    any(X(Fpos(i, Filt),:),1);
                
                % For an OR with positive output to become negative the
                % least needed change is the maximum of all positve inputs
                sF = RowValues(i,:);
                RowMargins(i, sF) = max(X(Fpos(i, Filt),sF) .* Xmargins(Fpos(i, Filt),sF),[],1);
                
                % For an OR with negative output to become positive the
                % least needed change is the minimum of all inputs (because they are all negative)
                sF = ~RowValues(i,:);
                RowMargins(i,sF) = min(Xmargins(Fpos(i, Filt),sF),[],1);
            else
                nFilt = Fneg(i,:)>0;
                if any(nFilt)
                    %% handle the NOT part
                    RowValues(i, :) =       ~X(Fneg(i, nFilt),:);
                    
                    % For an NOT to flip output is straightforward. Here I
                    % assume there is not more than one NOT in a branch.
                    RowMargins(i,:) = Xmargins(Fneg(i, nFilt),:);
                end
            end
        end
        Result(K,:) = all(RowValues,1);
        
        % Margin for samples that are classified as posive, is the minimum
        % fliping margin of all branches (because they are all positive)
        sF = Result(K,:);
        Margin(K,sF) = min(RowMargins(:,sF),[],1);
        
        % Margin for samples that are classified as negative, is the
        % maximum fliping margin of all negative branches
        RowMargins(RowValues) = 0; % To exclude positive branches
        sF = ~Result(K,:);
        Margin(K,sF) = max(RowMargins(:,sF),[],1);
        
    end
end

end