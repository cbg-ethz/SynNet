function Problem = Fit_Continuous_Classifier(Problem)
global Consts
try
    Plotting = Consts.Plotting;
    fprintf('Optimization started.\n');
    addpath('Auxilary_Functions/');
    TIC = tic;
    %% General  parameters
    MaxHitmenRatio        = Consts.Learning_MaxPoolRatio;
    NTopHitsLim           = Consts.Learning_MaxHitsToReport;
    MaxRounds             = Consts.Learning_MaxOptimizationRounds;
    Conv_THR1 = Consts.Learning_Convergence_Thr_1;  % after this many times of having same best performance program stops
    Conv_THR2 = Consts.Learning_Convergence_Thr_2;  % if AUC==100 after this many times of having same best performance program stops
    MinComplexity         = .5; % This is the percentage of unique unredundant circuits in the top of the pool. If it goes less than this I refine the whole top half of the pool
    %% Extract Params from the input structure
    MaxAnd    = Consts.MaxAnd;% Maximum number of AND inputs
    MaxOr     = Consts.MaxOr; % Maximum number of OR  inputs
    
    GenesData = Problem.Data.Values; % Data matrix, with samples on rows, and genes on columns
    Annotation = Problem.Data.Annots; % Sample annotation, Cancer=1, healthy=0
    Labels    = Problem.Data.GenIDs; % Gene Names
    SLabels   = Problem.Data.Labels; % Sample Names
    Nsamples  = size(GenesData,2);
    Ngene     = size(GenesData,1);
    Elites2Keep = max(10, NTopHitsLim * 3);
    MaxPoolSize = max([MaxHitmenRatio * Ngene, Elites2Keep*(1/MinComplexity)*10, 100]);
    MaxPoolSize = ceil(MaxPoolSize/2)*2; % This is the largest the pool of classifier circuits is going to get
    MaxEliteNum = ceil(.1*MaxPoolSize/2)*2;
    
    PosN =  sum( Annotation);
    NegN =  sum(~Annotation);
    PosI = find( Annotation);
    NegI = find(~Annotation);
    
    Labels_short = strrep(lower(Labels)        , 'hsa-', '');
    %% Quantize expressions and build a Recruit Index for each sample
    % Quanitization is not really needed for this it's just done for
    % the end report when we want to produce also binary statistics.
    [Problem.Data.BValues, ~ , Problem.Data.BMargin, Problem.Data.Threshold] = Quantize_Expresison(GenesData);
    Consts.BMargin_Max        = max(Problem.Data.BMargin(:));
    Consts.BMargin_Coeficient = 1./ (length(Annotation) * Consts.BMargin_Max);% To make sure it's always smaller that the smallest posisble change due to change in classifier accuracy
    
    Fresh_Recruits = BuildRecruitIndex_C(Problem);
    %% Initialization
    CPool.Function = [];% Initialize an empty Pool of functions
    tmpNprn=0;MaxStuckLog=0;
    fprintf(['Round\tTop performance\n']);
    Round = 0; NotConverged = true; AvgTopPerf= -1;Stuck =0;
    if Plotting
        %% Initialize optmization internal Variables plot
        Fig = figure('Color',[1 1 1]);
        subplot(2,6,[1:5 7:11])
        title(Problem.CancerLbl)
        hold on; xlabel('Round'); ylabel('Percentage'); ylim([-1 101]); grid on; box on
        SuccRatioold = 1;PoolComplexityold=100;AvgTopPerfold=zeros(NTopHitsLim,1);bestAUCold = zeros(NTopHitsLim,1);bestMarginold = zeros(NTopHitsLim,1);
        plot([Round], 10.^AvgTopPerfold(1)  , '-d', 'Color',[.0 .5 .15],'MarkerFaceColor' , [.0 .5 .15], 'linewidth', 1,'DisplayName', 'Top AUCs'     , 'markersize', 4);
        plot(Round, 0    , '-s', 'Color',[1 .5 .1],'MarkerFaceColor' , [1 .5 .1], 'linewidth', .5,'DisplayName', 'Top Avg. Margin'     ,  'markersize', 4);
        plot([Round], [PoolComplexityold]   , '-s', 'Color',[.5 .5 .50],'MarkerFaceColor' , [.5 .5 .50], 'linewidth', 1,'DisplayName', 'Function Pool Complexity'   , 'markersize', 4);
        plot([Round], [SuccRatioold ]*100   , '-s', 'Color',[.0 .25 .55],'MarkerFaceColor' , [.0 .25 .55], 'linewidth', 1,'DisplayName', 'NewGene Success Rate'       , 'markersize', 4);
        
        Lh = legend('show', 'location','NorthWest');
        set(Lh, 'Interpreter','tex')
        set(gcf, 'position', [200 200 1200 500]);
        set(gca,'Color',[.0 .0 0]);
        set(gcf, 'InvertHardCopy', 'off');
        
        subplot(2,6,6)
        axis off
        xlim([-1, 1])
        ylim([-1, 1])
        BestCirc_Y = 1;
        if isfield(Problem, 'Synthesis_Function')
            [~, ~, Continuous_Stats] = Get_Absolute_Performance_C(Problem.Synthesis_Function, GenesData, Annotation);
            FunString = SPrint_Function(Problem.Synthesis_Function, Labels_short, [Continuous_Stats.AUC, Continuous_Stats.Margin]);
            FunString = sprintf('Synthesis Function:\n%s\n', FunString);
            text(-1,1,FunString, 'horizontalalignment', 'left', 'verticalalignment', 'top', 'fontsize', 10);
            BestCirc_Y = 0;
        end
        pT = text(0,0,'');
        hold on
    end
    %% Iterative optimization
    while NotConverged && MaxRounds > Round
        
        CleanUpRound = false;
        Round = Round + 1;
        %% Fill the Classifier hitpool if empty
        if isempty(CPool.Function)
            CPool.Function         = Fresh_Recruits.Function;
            tmpPS = size(CPool.Function,1);
            CPool.IsARookie        = ones(tmpPS,1);
            CPool.IsFormatted      = true(tmpPS,1);
            CPool.IsNew            = true(tmpPS,1); clear tmpPS
        end
        
        %% Evaluate pool: Hunt for cancer
        CPool.WPerformance(CPool.IsNew,1) = nan; % Just to initiate it
        [CPool.APerformance(CPool.IsNew,1), CPool.Log10_Function_Output(CPool.IsNew,:)] = Get_Absolute_Performance_C(CPool.Function(CPool.IsNew,:,:), GenesData, Annotation);
        %% Mark Elite Classifiers, and sort the rest by weighted score
        [Elites, nElites] = Catch_Elites(CPool.APerformance, Elites2Keep);
        uEliteI = Remove_ConsecutiveRepeats(CPool.Function(Elites,:,:), CPool.APerformance(Elites));
        tmpBestHitMan = CPool.Function(Elites(uEliteI),:,:);
        BestHitMan      = tmpBestHitMan(1:NTopHitsLim,:,:);
        PoolComplexity  = size(tmpBestHitMan,1) / nElites;
        nElites = size(tmpBestHitMan,1);
        CutFrom = max(nElites, min(MaxPoolSize/2, size(CPool.Function,1))); % Discard up to 50% of the current pool based on performance. Keep at least 3 times more than NTopHitsLim
        
        CPool.WPerformance = Get_Weighted_Performance(CPool.Log10_Function_Output, Annotation, Elites(uEliteI));
        
        [~, WperfI]  = sort(CPool.WPerformance, 'descend');
        CPool = Pej_Struct_RowSelect(CPool, WperfI); % Sort by performance
        CPool.IsNew(:) = false;
        
        Rp = (1 + sum(CPool.IsARookie(1:CutFrom)==1)) / (2+sum(CPool.IsARookie==1));
        Rr = (1 + sum(CPool.IsARookie(1:CutFrom)==0))/  (2+sum(CPool.IsARookie==0));
        SuccRatio = (Rp)/(Rp+Rr);
        
        if (PoolComplexity < MinComplexity) || (nElites >= MaxEliteNum)
            %% Clear the pool from repeats, if there are too many equivalent classifiers in the hitpool
            iRfP = Remove_ConsecutiveRepeats(CPool.Function, CPool.APerformance);
            CPool = Pej_Struct_RowSelect(CPool, iRfP);
            CleanUpRound = true;
            CutFrom = min(CutFrom, size(CPool.Function,1));
            
        else
        end
        if 1%~CleanUpRound
            %% Reporting, and convergence check
            [BestHitManPerf, ~, Continuous_Stats] = Get_Absolute_Performance_C(BestHitMan, GenesData, Annotation);
            NewAvgTopPerf   = (BestHitManPerf);
            bestAUC = Continuous_Stats.AUC;
            bestMargin = Continuous_Stats.Margin;
            
            if all(AvgTopPerf == NewAvgTopPerf)
                % It's suck/converged
                Stuck = Stuck + 1;
            elseif any(AvgTopPerf > NewAvgTopPerf)
                % We lost one of our best solutions! This should not happen.
                fprintf('\nWARNING: Performance decreased!\n')
                tmpNprn = 0;
                Stuck = Stuck + 1;
            else
                % It Improved Again!
                Stuck = 0;
            end
            
            if bestAUC(1) <100
                % We have not found yet any correct order
                MaxStuckLog = max(Stuck,MaxStuckLog); % This reported. It's useful, in hindsight, in estimating the confidence of method based on the "Conv_THR1" value.
            end
            
            Secondary_Convergence = (bestAUC(NTopHitsLim)==100) && (Stuck >= Conv_THR2);
            NotConverged = (Stuck < Conv_THR1) && (~Secondary_Convergence);
            
            AvgTopPerf = NewAvgTopPerf;
            
            fprintf(repmat('\b', 1, tmpNprn));
            if Plotting
                subFig = subplot(2,6,12);
                [~, ~, ~] = Get_Continous_Margins(BestHitMan(1,:,:), GenesData, Annotation, subFig);
            else
                %                 [~, Continous_MarginW, Continous_MarginA] = Get_Continous_Margins(BestHitMan(1,:,:), GenesData, Annotation);
            end
            tmpNprn = fprintf('%d\t Current best stats: AUC %.2f\tMargin %.2f(log10)\tStuck %d(rounds)', Round, bestAUC(1), bestMargin(1), Stuck);
            
            if Plotting
                %% Plot optmization internal Variables
                figure(Fig);
                subplot(2,6,[1:5 7:11])
                Prx = [Round-1 Round];
                plot(Prx, [SuccRatioold SuccRatio]*100           , '-s', 'Color', [.0 .25 .55],'MarkerFaceColor' ,  [.0 .25 .55], 'linewidth', .5,'DisplayName', 'NewGene Success Rate'       , 'markersize', 2);
                plot(Prx, 100*[PoolComplexityold PoolComplexity] , '-s', 'Color',[.4 .4 .45],'MarkerFaceColor' , [.4 .4 .45], 'linewidth', .5,'DisplayName', 'Is formatted rate'          ,  'markersize', 2);
                if CleanUpRound
                    plot([Round Round], [PoolComplexity 1]*100, '--','Color',[.4 .4 .45],'MarkerFaceColor' , [.0 .0 .15],  'linewidth', 1 , 'markersize', 4);
                    plot(Round        ,    PoolComplexity *100, 's' ,'Color',[1   0   0],'MarkerFaceColor' , [1   0   0],  'linewidth', 1 , 'markersize', 4);
                    PoolComplexity = 1;
                end
                
                plot(Prx, 10.^[bestMarginold(1) bestMargin(1)], '-s' , 'Color',[1 .5 .1],'MarkerFaceColor' , [1 .5 .1], 'linewidth', .5,'DisplayName', 'Continuous Margin'     ,  'markersize', 4);
                plot(Prx,         [bestAUCold,    bestAUC] , ':' , 'Color',[.0 .5 .15],'MarkerFaceColor' , [.0 .5 .15], 'linewidth', .5,'DisplayName', 'Best CPool.Performance'     ,  'markersize', 4);
                plot(Prx,         [bestAUCold(1) bestAUC(1)]  , '-s', 'Color',[.0 .5 .15],'MarkerFaceColor' , [.0 .5 .15], 'linewidth', .5,'DisplayName', 'Best CPool.Performance'     ,  'markersize', 4);
                
                xlim([Round-75, Round+25]);
                drawnow();
                SuccRatioold        = SuccRatio;
                PoolComplexityold   = PoolComplexity;
                bestAUCold          = bestAUC;
                bestMarginold      = bestMargin;
                % Plot the current best function
                subplot(2,6,6)
                FunString = SPrint_Function(BestHitMan(1,:,:), Labels_short, [bestAUC(1), bestMargin(1)]);
                FunString = sprintf('\n\n\nPool size:\t%d,\tbest function so far:\n%s', size(CPool.Function,1), FunString);
                delete(pT);
                pT = text(-1, BestCirc_Y,FunString, 'horizontalalignment', 'left', 'verticalalignment', 'top', 'fontsize', 10);
            end
            
            %% Exit if converged
            if ~NotConverged
                break
            end
            
            %% Improve
            nNewHits = MaxPoolSize - CutFrom;
            EliteKids = round(nNewHits/2); NormalKids = nNewHits - EliteKids;
            Parents = [randi(nElites, EliteKids, 1); randi(CutFrom-nElites, NormalKids, 1)+nElites]; % Use the Elites and normal hits 50/50% to make the next generation
            
            NewHits = CPool.Function(Parents, :,:);NewHitsRookies = ones(nNewHits, 1);
            
            if mod(Round,2)
                %% Improve Negative samples
                S2improve = NegI(randi(NegN,nNewHits, 1)); % I choose this sample to be improved.
            else
                %% Improve positive samples
                S2improve = PosI(randi(PosN,nNewHits, 1)); % I choose this sample to be improved.
            end
            for h = nNewHits:-1:1
                if ~CleanUpRound && (rand > .5)%SuccRatio)
                    % Recombination
                    RecombPart  = randi(CutFrom);
                    FullBranches = find(any(CPool.Function(RecombPart,:,:),3));
                    numBranches = length(FullBranches);
                    if numBranches ==0
                        BranchI = 1; % Basically this will put zero somewhere.
                    else
                        BranchI = FullBranches(randi(numBranches));
                    end
                    TragetBranchI = randi(MaxAnd);
                    TargetBranch  = NewHits(h,TragetBranchI,:);
                    if any(TargetBranch<0)
                        NewHits(h,TragetBranchI,:) = CPool.Function(RecombPart,BranchI,:);
                    else
                        TwoBranches = [TargetBranch, CPool.Function(RecombPart,BranchI,:)];
                        NewHits(h,TragetBranchI,:) = TwoBranches(randperm(MaxOr*2,MaxOr));
                    end
                    NewHitsRookies(h) = 0;
                else
                    % use a new
                    sHandler = find(Fresh_Recruits.SampleWise_WeightCumSum(:,S2improve(h)) >= rand,1, 'first'); % Choose one of the genes that can do sample "s".
                    if Fresh_Recruits.IsaNOT(sHandler)
                        % it's a NOT it can't share a row
                        NewHits(h,randi(MaxAnd),:) = [Fresh_Recruits.GeneID(sHandler), zeros(1, MaxOr-1)]; % Overwrite a row
                    else
                        % It's an OR can be combined in the lines
                        NewHits(h,randi(MaxAnd),randi(MaxOr)) = Fresh_Recruits.GeneID(sHandler); % Throw it somehwere (if it conflicts with a NOT the branch will be removed anyway)
                    end
                    
                end
            end
            
            [~, NewHits] = Ispractical(NewHits,'Prune');
            NewHits = FixFunctionFormat(NewHits);
            if size(NewHits,1)+CutFrom < size(CPool.Function,1)
                CPool = Pej_Struct_RowDel(CPool, (CutFrom+size(NewHits,1)+1):size(CPool.Function,1));
            end
            CPool.IsNew  = [false(CutFrom, 1); true(nNewHits, 1)];
            
            CPool.Function(CPool.IsNew,:,:) = NewHits;
            CPool.IsARookie(~CPool.IsNew) = nan;
            CPool.IsARookie(CPool.IsNew) = NewHitsRookies;
            CPool.IsFormatted(CPool.IsNew)  = true;
            
        end
        
    end
    fprintf('\nBest Circuit Found:\n\n%s\n', SPrint_Function(BestHitMan, Labels, [bestAUC, bestMargin]));
    
    toc(TIC)
    %% Report results
    Problem.Results.BestCircuits      = BestHitMan;
    Problem.Results.BestCircuitsPerfs = BestHitManPerf;
    Problem.Results.ExitFlags.HitMaxRounds  = MaxRounds <= Round;
    Problem.Results.ExitFlags.WasStuck      = MaxStuckLog;
    Problem.Results.ExitFlags.ElapsedTime   = toc(TIC);
    
    Problem.Consts = Consts;
    if Plotting
        Pej_SavePlot(Fig, [Consts.OutFolder '/FitPlot_' Problem.CancerLbl]);
    end
catch Err
    disp('Error!')
    Err.getReport
end
end