function Plot_binary_output(Function, Data, Annotation, FigurePath)
nAnnotation = ~Annotation;
[Function_Output, Output_Margins]  = Evaluate_Function_B_Margin(Function, Data.BValues, Data.BMargin);
isNotCorrect = Function_Output ~= Annotation;

Output_Margins = min(Output_Margins, 2)/2; % So I don't care about margins greater than 10^2, it's good enough! and I normalize them to be less than .5
% Output_Margins(:,nAnnotation) = - Output_Margins(:,nAnnotation);
Output_Margins(:,~Function_Output) = - Output_Margins(:,~Function_Output);
% Output_Margins(isNotCorrect& Annotation) = -.5;
% Output_Margins(isNotCorrect&nAnnotation) = +.5;

if ~ischar(FigurePath)
    subplot(FigurePath);
    % Assume it's a figure
    Sp = sum(Annotation);Sn = sum(nAnnotation);
    Xp = linspace(0,.8,Sp); Xp = Xp - mean(Xp) + .55;
    Xn = linspace(0,.8,Sn); Xn = Xn - mean(Xn) - .55;
    RND(Annotation)  = Xp;
    RND(nAnnotation) = Xn;
    xlim([-1 1]);
    hold off
    plot(xlim, [0 0], '-k')
    box on
    hold on
    plot(RND(Annotation), Output_Margins(Annotation)...
        ,'s', 'markersize', 5, 'MarkerFaceColor' , [.5 .05 .1], 'displayname', 'Positive', 'markeredgecolor', [.3 .05 .1])
    plot(RND(nAnnotation), Output_Margins(nAnnotation)...
        , 's', 'markersize', 5, 'MarkerFaceColor' , [.1 .05 .6], 'displayname', 'Negative','markeredgecolor', [.1 .05 .3])
    hold off
    set(gca, 'Xtick', [-.5 .5])
    set(gca, 'XTIckLabel', {'Pos' 'Neg'})
    ylim([-1.05 1.05])
    set(gca, 'Ytick', [-1 -.5 -.1 .1 .5 1])
    set(gca, 'YtickLabel', {'>2' '1' 'Off', 'On' '1' '>2'})
    ylabel('Margin (log_{10})')
else
    % Assume it's a path for figure
    Fig = figure;
    
    RND = rand(size(Annotation))*.4 -.2;
    xlim([-.5 1.5])
    plot(xlim, [0 0], '-k')
    hold on
    plot(Annotation(Annotation)+RND(Annotation), Output_Margins(Annotation)...
        ,'s', 'markersize', 6, 'MarkerFaceColor' , [.5 .05 .1], 'displayname', 'Positive', 'markeredgecolor', [.3 .05 .1])
    plot(Annotation(nAnnotation)+RND(nAnnotation), Output_Margins(nAnnotation)...
        , 's', 'markersize', 6, 'MarkerFaceColor' , [.1 .05 .6], 'displayname', 'Negative','markeredgecolor', [.1 .05 .3])
    %legend show
    set(gca, 'Xtick', [0 1])
    set(gca, 'XTIckLabel', {'Healthy' 'Cancer'})
    box on
    ylabel('Circuit output (mol/cell)')
    set(gcf, 'position', [200 200 330 330]);
    ylim([-1.05 1.05])
    set(gca, 'Ytick', [-1 -.5 -.1 .1 .5 1])
    set(gca, 'YtickLabel', {'>2' '1' 'Off', 'On' '1' '>2'})
    ylabel('Margin (log_{10})')
    
    
    tmpMp  = mean(Function_Output( Annotation));
    tmpMn  = mean(Function_Output(nAnnotation));
    tmpM   = mean([tmpMp; tmpMn]);
    Perf = (tmpMp - tmpMn + 1)/2; % (informedness + 1 )/2
    text(.5, .6, sprintf('%.1f%%', Perf*100), 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'BackgroundColor', [1 1 1], 'fontsize', 14);
    errorbar([.5 5], [0 5], [.5, 0], '.', 'linewidth', 1.2, 'color', [.1 .1 .1])
    ylabel('Margin (log_{10})')
    xlim([-.5 1.5])
    Pej_SavePlot(Fig, FigurePath)
end
