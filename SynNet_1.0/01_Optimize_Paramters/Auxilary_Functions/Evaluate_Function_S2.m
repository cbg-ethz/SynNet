% This file evaluates the continous for of a boolean function over values in "X",

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
%   IsVagueX(leftover from the binary version): THIS SHOULD BE LEFT OUT OR SET TO ALL ZEROS, I just keep the parameter
%   to keep the consistency with the binary version. This is a binary matrix, the same size as "X",
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

function [Result IsVague RowValues IsVagueRow] = Evaluate_Function_S2(FunctionArray ,Xorg, StepPoint)
t = log10(Xorg)-log10(StepPoint);
X = 1./(1+exp(-t*20));
X(X>(1-1E-6))=1;
X(X<(  1E-6))=0;
clear Xorg
IsVagueX = false(size(X)); % If there's no info on the vagueness of dta, assume nothing's vague.
if any(IsVagueX); warning('VAGUE Value in the data! Theres is no such thing in the continous evaluation!'); end
if length(size(FunctionArray))==2
    %% It's a Single function, format it into an 3D array of length one.
    tmpFnc = FunctionArray;
    clear FunctionArray
    FunctionArray(1,:,:) = tmpFnc;
    clear tmpFnc
end
MaxAnd = size(FunctionArray,2); % Number of AND inputs in the function family
N      = size(X,2);             % Number of Samples in the data.

for K = size(FunctionArray,1):-1:1
    RowValues  = zeros(MaxAnd, N);
    FullLine   = true(MaxAnd,1);
    Function(:,:) = (FunctionArray(K,:,:));
       %% handle the OR part
    for i = MaxAnd:-1:1
        Line = Function(i,:);
        if any(Line<0)
            %% NOT
            RowValues(i,:) = 1-X(-sum(Line),:); % I ASSUME THERE's ONLY ONE NOT IN A LINE
        elseif  any(Line>0)
            %% OR
            Fpos = Line(Line>0);
            RowValues(i,:) = max(X(Fpos,:),[],1);
        else
            % its empty
            FullLine(i)=false;
        end
    end
    RowValues(~FullLine, :)=[];
    Result(K,:) = min(RowValues,[],1);
end


if nargout >1; IsVague      = false;  end
if nargout >2; RowValues    = 0;      end
if nargout >3; IsVagueRow   = false;  end

end
