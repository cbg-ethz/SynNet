function [MarginM, MarginW, MarginA, RoC_AUC] = Get_Continous_Margins(Function, GenesData, TrueAnnoation, FigurePath)
AnnoationCnt = Evaluate_Function_C(Function , GenesData);
MarginW = Discrimination_Margin_Worst(AnnoationCnt, TrueAnnoation);
MarginM = Discrimination_Margin_Median(AnnoationCnt, TrueAnnoation);
MarginA = Discrimination_Margin_Mean(AnnoationCnt, TrueAnnoation);
if nargout >3
    %     [~, ~, ~, RoC_AUC] = perfcurve(TrueAnnoation, AnnoationCnt, true);
    RoC_AUC =    fastAUC(double(TrueAnnoation), AnnoationCnt,  1);
end

if nargin > 3
    if ~ischar(FigurePath)
        subplot(FigurePath);
        % Assume it's a figure
        Sp = sum(TrueAnnoation);Sn = sum(~TrueAnnoation);
        Xp = linspace(0,.8,Sp); Xp = Xp - mean(Xp) + .55;
        Xn = linspace(0,.8,Sn); Xn = Xn - mean(Xn) - .55;
        RND(TrueAnnoation)  = Xp;
        RND(~TrueAnnoation) = Xn; 
        semilogy(RND(TrueAnnoation), AnnoationCnt(TrueAnnoation)...
            ,'s', 'markersize', 5, 'MarkerFaceColor' , [.5 .05 .1], 'displayname', 'Positive', 'markeredgecolor', [.3 .05 .1])
        hold on
        semilogy(RND(~TrueAnnoation), AnnoationCnt(~TrueAnnoation)...
            , 's', 'markersize', 5, 'MarkerFaceColor' , [.1 .05 .6], 'displayname', 'Negative','markeredgecolor', [.1 .05 .3])
        hold off
        xlim([-1 1]);
        set(gca, 'Xtick', [-.5 .5])
        set(gca, 'XTIckLabel', {'H' 'C'})
    else
        % Assume it's a path for figure
        Fig = figure;
        RND = rand(size(TrueAnnoation))*.4 -.2;
        semilogy(TrueAnnoation(TrueAnnoation)+RND(TrueAnnoation), AnnoationCnt(TrueAnnoation)...
            ,'s', 'markersize', 6, 'MarkerFaceColor' , [.5 .05 .1], 'displayname', 'Positive', 'markeredgecolor', [.3 .05 .1])
        hold on
        semilogy(TrueAnnoation(~TrueAnnoation)+RND(~TrueAnnoation), AnnoationCnt(~TrueAnnoation)...
            , 's', 'markersize', 6, 'MarkerFaceColor' , [.1 .05 .6], 'displayname', 'Negative','markeredgecolor', [.1 .05 .3])
        %legend show
        xlim([-.5 1.5])
        set(gca, 'Xtick', [0 1])
        set(gca, 'XTIckLabel', {'Pos' 'Neg'})
        box on
        ylabel('Circuit output (mol/cell)')
        set(gcf, 'position', [200 200 330 330]);
        grid on
        if max(ylim)/min(ylim) < 10
            yL = log10(ylim);
            yL(1)= floor(yL(1));
            yL(2)=  ceil(yL(2));
            ylim(10.^yL)
        end
        
        tmpMp  = geomean(AnnoationCnt( TrueAnnoation));
        tmpMn  = geomean(AnnoationCnt(~TrueAnnoation));
        tmpM   = geomean([tmpMp; tmpMn]);
        expFc  = 10.^MarginA;
        text(.5, tmpMp*1.1, sprintf('%.1f', expFc), 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'BackgroundColor', [1 1 1], 'fontsize', 14);
        errorbar([.5 5], [tmpM 5], [tmpM - tmpMn, 0], [tmpMp - tmpM, 0], '.', 'linewidth', 1.2, 'color', [.1 .1 .1])
        xlim([-.5 1.5])
        Pej_SavePlot(Fig, FigurePath)
    end
end