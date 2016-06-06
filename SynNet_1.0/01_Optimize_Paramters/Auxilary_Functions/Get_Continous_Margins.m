function [MarginM, MarginW] = Get_Continous_Margins(Function, GenesData, TrueAnnoation, FigurePath)
AnnoationCnt = Evaluate_FunctionCnt(Function , GenesData);
MarginW = Discrimination_Margin_Worst(AnnoationCnt, TrueAnnoation);
MarginM = Discrimination_Margin_Median(AnnoationCnt, TrueAnnoation);

if nargin > 3
   RND = rand(size(TrueAnnoation))*.4 -.2;
    Fig = figure;
   semilogy(TrueAnnoation(TrueAnnoation)+RND(TrueAnnoation), AnnoationCnt(TrueAnnoation)...
       ,'s', 'markersize', 6, 'MarkerFaceColor' , [.5 .05 .1], 'displayname', 'Cancer', 'markeredgecolor', [.3 .05 .1])
   hold on
   semilogy(TrueAnnoation(~TrueAnnoation)+RND(~TrueAnnoation), AnnoationCnt(~TrueAnnoation)...
       , 's', 'markersize', 6, 'MarkerFaceColor' , [.1 .05 .6], 'displayname', 'Healthy','markeredgecolor', [.1 .05 .3])
   %legend show
   xlim([-.5 1.5])
   set(gca, 'Xtick', [0 1])
   set(gca, 'XTIckLabel', {'Healthy' 'Cancer'})
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
   Pej_SavePlot(Fig, FigurePath)
end