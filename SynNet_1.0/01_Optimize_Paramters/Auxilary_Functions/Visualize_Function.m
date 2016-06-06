function Visualize_Function(FunctionArray, Problem)
OutDir  = ['Optimization_Figures_' num2str(round(Problem.IdealFunction.stepAt)) '/'];
mkdir(OutDir);
global Functioninstance Domain oBP

Step_High = Problem.IdealFunction.Step_High;
Step_Low  = Problem.IdealFunction.Step_Low;
InitialParams = Problem.Continuous_Function.InitialParams;
OPTtheta = Problem.OPTtheta;
eOPTtheta = 10.^OPTtheta;
eInitialParams = 10.^InitialParams;


for f = 1: size(FunctionArray,1)
    %% Distribution Plot
    Functioninstance = zeros(size(FunctionArray,2), size(FunctionArray,3));
    Functioninstance(:,:) = FunctionArray(f, :,:);
    Fig = figure;
    set(gcf, 'position', [100 100 300 300]);
    xlim([-1 5])
    Xs = DrawSamples(10000);
    Initial_Output  = (Evaluate_FunctionCnt_Minimal(Xs, eInitialParams))';
    Ideal_Output    = IdealFun(Xs)>((Step_High -Step_Low)/2);
    Optimized_Output= (Evaluate_FunctionCnt_Minimal(Xs, eOPTtheta))';
    set(gca, 'YScale', 'log');
    
    % Plot Initial distributions
    tmpPos = Initial_Output( Ideal_Output);
    tmpNeg = Initial_Output(~Ideal_Output);
    tmpMp  = geomean(tmpPos);
    tmpMn  = geomean(tmpNeg);
    tmpM   = geomean([tmpMp; tmpMn]);
    expFc  = log10(tmpMp)- log10(tmpMn);
    text(3.5, tmpMp*1.1, sprintf('%.2f', expFc), 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'BackgroundColor', [1 1 1])
    Pej_Plot_Dist_Customized([2.8 2.8], repmat(tmpPos,1,2),  [.7 .05 .1],[-1; 1]);%  ,...
    Pej_Plot_Dist_Customized([4.2 4.2], repmat(tmpNeg,1,2),  [.7 .05 .1],[-1; 1]);%  ,...
    errorbar(3.5, tmpM, tmpM - tmpMn, tmpMp - tmpM, 'linewidth', 1.2, 'color', [.1 .1 .1]);
    % Plot Optimized distributions
    tmpPos = Optimized_Output( Ideal_Output);
    tmpNeg = Optimized_Output(~Ideal_Output);
    tmpMp  = geomean(tmpPos);
    tmpMn  = geomean(tmpNeg);
    tmpM   = geomean([tmpMp; tmpMn]);
    expFc  = log10(tmpMp)- log10(tmpMn);
    text(0.5, tmpMp*1.1, sprintf('%.2f', expFc), 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'BackgroundColor', [1 1 1])
    errorbar([0.5 3.5], [tmpM nan], [tmpM - tmpMn nan], [tmpMp - tmpM nan], 'linewidth', 1.2, 'color', [.1 .1 .1]);
    Pej_Plot_Dist_Customized([-.2 -.2], repmat(tmpPos,1,2),  [.1 .05 .5],[-1; 1]);%  ,...
    Pej_Plot_Dist_Customized([1.2 1.2], repmat(tmpNeg,1,2),  [.1 .05 .5],[-1; 1]);%  ,...
    
    
    
    xlim([-1 5])
    box on
    plot([2 2], ylim, 'k', 'linewidth', 1)
    ylabel('Output marker distribution(mol/cell)')
    set(gca, 'xtick', [0 1 3 4])
    set(gca, 'xtickLabel', {'Cancer' 'Healthy' 'Cancer' 'Healthy'})
    grid on
    CName= SPrint_Function(Functioninstance);
    CName(CName==10)=[];% Remove line breaks
    title(CName)
    Pej_SavePlot(Fig, [OutDir 'Function_' int2str(f)]);
    
    if max(abs(Functioninstance(:)))==1
        %% Operation Curve for 1D functions
        Xs  = logspace(Domain(1,1), Domain(1,2), 50);
        Ix  = Evaluate_Function_Minimal(Xs>oBP);
        Fxi = Evaluate_FunctionCnt_Minimal(Xs, eInitialParams);
        Fxo = Evaluate_FunctionCnt_Minimal(Xs, eOPTtheta);
        
        Fig = figure;
        set(gcf, 'position', [100 100 400 175]);
        loglog(Xs, Fxi, ...
            's-', 'markersize',  3, 'markerfacecolor',[.6 .05 .1], 'color',[.6 .05 .1], 'Displayname', 'Original Function');hold on;
        semilogx(Xs, Ix * max(Fxo) + min(Fxo)  ,...
            's-', 'markersize',  3, 'markerfacecolor',[.01 .5 .1],'color',[.1 .5 .1], 'Displayname', 'Ideal Function');
        semilogx(Xs, Fxo,...
            's-', 'markersize',  3, 'markerfacecolor',[.1 .05 .5],'color',[.1 .05 .5],'Displayname', 'Optimized Function');
        legend('show', 'location', 'Best')
        xlim(10.^(Domain(1,:)))
        ylabel('Output marker (mol/cell)')
        xlabel('Input Expression (mol/cell)')
        title([CName])
        %         grid on
        Pej_SavePlot(Fig, [OutDir 'Function_' int2str(f) 'Curv']);
        
    end
    if max(abs(Functioninstance(:)))==2
        %% Surface plot for 2D functions
        [X, Y] = meshgrid(logspace(Domain(1,1), Domain(1,2), 50), logspace(Domain(2,1), Domain(2,2), 50));
        XY = [X(:)'; Y(:)'];
        Zi = reshape(Evaluate_FunctionCnt_Minimal(XY, eInitialParams), size(X));
        Zi = Transform(Zi);
        Zo = reshape(Evaluate_FunctionCnt_Minimal(XY, eOPTtheta), size(X));
        Zo = Transform(Zo);
        CRange = [min([min(Zi(:)), min(Zo(:))]), max([max(Zi(:)), max(Zo(:))])];
        Fig = figure;
        set(gcf, 'position', [100 100 700 300]);
        %         xlabel('G1 expression')
        %         ylabel('G2 expression')
        subplot(1, 7,[1 3])
        hold on
        set(gca, 'Yscale', 'log');
        set(gca, 'Xscale', 'log')
        H = pcolor(X, Y, Zi);
        set(H,'linewidth',.1);
        colormap jet
        
        caxis(CRange);
        axis square
        title([CName ' (original)'])
        set(gca, 'position', [ 0.1300    0.1100    0.3167    0.8150])
        
        
        subplot(1, 7,[4 7])
        hold on
        set(gca, 'Yscale', 'log');
        set(gca, 'Xscale', 'log')
        H = pcolor(X, Y, Zo);
        set(H,'linewidth',.1);
        colormap jet
        caxis(CRange);
        axis square
        title([CName ' (optimized)'])
        set(gca, 'ytick', [])
        %         xlabel('G1 expression')
        
        %         subplot(1, 7,[7])
        %         axis off
        %         caxis(CRange);
        colorbar
        set(gca, 'position', [ 0.4739    0.1100    0.3167    0.8150])
        Pej_SavePlot(Fig, [OutDir 'Function_' int2str(f) 'Surf']);
        
    end
end

Fig = figure;
set(gcf, 'position', [100 100 300 300]);
Params2report = [OPTtheta, InitialParams];
hold on
bar(Params2report);
errorbar(1:length(Params2report), InitialParams, InitialParams - Problem.lb, Problem.ub - InitialParams , '+', 'linewidth', 1.5, 'color',  [.01 .75 .1])
text(1:length(Params2report), zeros(size(InitialParams)), {'C_{1} '; 'C_{2} '; 'T_{max} '; 'FF4_{max} '; 'Out_{max} '}, 'horizontalalignment', 'right', 'fontsize', 8,'rotation', 90)
ylim([-max(ylim)*.2 max(ylim)])
set(gca, 'XTick', []);
ytmp = get(gca, 'YTick');
set(gca, 'YTick', ytmp(ytmp>=0));

ylabel('mol/cell (log_{10})')
grid on
box on

Pej_SavePlot(Fig,  [OutDir 'Parameters'])


end

function Xs = DrawSamples(N)
%Make <=N samples from the background distribution of the data
global PDname Pa Pb Pc Domain
Filt = false(1,N);
for d = size(Domain,1):-1:1
    if isnan(Pb)
        Xs(d,:) = random(PDname{d}, Pa(d),1,N);
    elseif isnan(Pc)
        Xs(d,:) = random(PDname{d}, Pa(d), Pb(d), 1,N);
    else
        Xs(d,:) = random(PDname{d}, Pa(d), Pb(d), Pc(d), 1,N);
    end
    Filt = Xs(d,:)>Domain(d,2) | Xs(d,:)<Domain(d,1) | Filt;
end
Xs(:,Filt)=[];
Xs = 10.^(Xs);
end


function Z = Transform(Z)
Z = log10(Z);
Z = Z- min(Z(:));
% Z = Z-max(Z(:))/2;
end