% Warning: This Code is written assuming that the function would be small
% (~10 genes), it   uses exponential amount of memory and CPU!

function Pruned_Function = Prune_Circuit(Function, Data)
global Consts

if length(size(Function))==2
    %% It's a Single function, format it into an 3D array of length one, This will speed up the bootstrap
    tmpFnc = Function;
    clear FunctionArray
    Function = zeros(1,Consts.MaxAnd, Consts.MaxOr);
    Function(1,:,:) = tmpFnc;
    clear tmpFnc
end

if size(Function,1)>1
    % prune each of them separately
    for i = size(Function,1):-1:1
        Pruned_Function(i,:,:) = Prune_Circuit(Function(i,:,:), Data);
    end
else
    fprintf('Pruning ...\t')
    nBoot = 1000;
    A = Data.Annots;
    L = length(A);
    
    % I should also discard unnnecessary parts of the data!
    
    
    subF = Get_subFuncitons(Function);
    subF = Refine_HitPool(subF); % To remove logically redundant functions
    if Consts.AnalysisMode == 'B'
        D.BValues = Data.BValues;
        D.BMargin = Data.BMargin;
        for b = nBoot:-1:1
            BsI = randi(L, L,1);
            if all(A(BsI)) || ~any(A(BsI))
                % This case happens with the probability of the majority classs
                % to power of sample number
                bAUC(b) = .5;
                bMrg(b) = 0;
            else
                bD.BValues =  D.BValues(:,BsI);
                bD.BMargin =  D.BMargin(:,BsI);
                [~, ~, Binary_Stats] = Get_Absolute_Performance_B(Function, bD, A(BsI));
                bAUC(b) = Binary_Stats.RawPerformance;
                bMrg(b) = Binary_Stats.Margin;
            end
        end
        
        
        [~, ~, subF_Stats] = Get_Absolute_Performance_B(subF, D, A);
        P(:,1) = subF_Stats.RawPerformance;
        P(:,2) = subF_Stats.Margin;
    else
        D = Data.Values;
        for b = nBoot:-1:1
            BsI = randi(L, L,1);
            if all(A(BsI)) || ~any(A(BsI))
                % This case happens with the probability of the majority classs
                % to power of sample number
                bAUC(b) = .5;
                bMrg(b) = 0;
            else
                [~, ~, Continuous_Stats] = Get_Absolute_Performance_C(Function, D(:,BsI), A(BsI));
                bAUC(b) = Continuous_Stats.AUC;
                bMrg(b) = Continuous_Stats.Margin;
            end
        end
        
        
        [~, ~, subF_Stats] = Get_Absolute_Performance_C(subF, D, A);
        P(:,1) = subF_Stats.AUC;
        P(:,2) = subF_Stats.Margin;
    end
    
    madAUC = mad(bAUC,1);                % Median Abs. dev
    aCrit  = median(bAUC) - madAUC;      % Critical value for AUC, median - one MAD
    if madAUC == 0
        mCrit  = median(bMrg) - mad(bMrg,1); % Critical value for AUC, median - one MAD
    else
        % If there is a difference in the AUC of the models, we don't care about the margin!
        mCrit = -Inf;
    end
    
    % Remove all that are out of one STD
    Fa = P(:,1) < aCrit;
    Fm = P(:,2) < mCrit;
    F  = Fa | Fm;
    if all(F)
        % This can happen when the median of the bootstrap samples happens
        % to be just a little bit better than the actual network
        [~, tmpI_] = sortrows(P, [-1 -2]);
        F(tmpI_)= false;
    end
    P(F,:) = [];
    subF(F,:,:) = [];
    
    % remove those with higher complexitiy that the simplest one
    subF_nG = sum(sum(subF~=0,3),2);
    Fc  =  subF_nG >min(subF_nG);
    P(Fc,:) =[];
    subF(Fc,:,:) = [];
    
    % choose the best from what's left
    Fa  =  P(:,1) < max(P(:,1));
    P(Fa,:) =[];
    subF(Fa,:,:) = [];
    
    % choose the best from what's left
    Fm  =  P(:,2) < max(P(:,2));
    P(Fm,:) =[];
    subF(Fm,:,:) = [];
    
    Pruned_Function = subF(1, :, :);
    if size(subF,1)>1
        warning('there are more that one equally good pruned functions:\n %s', SPrint_Function(subF))
    end
    
    fprintf('done!\n')
end

end