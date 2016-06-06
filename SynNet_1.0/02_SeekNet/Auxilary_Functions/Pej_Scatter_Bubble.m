function Pej_Scatter_Bubble(plotX,plotY,plotZ, ColorFix)
PN = size(plotX,1);
if nargin < 4
    ColorFix = [];
    ColorBox = ones(PN,1)*[.1 .1 .7] + [rand(PN,2)*.2-.1 zeros(PN,1)];
    ColorBox = max(ColorBox,0);ColorBox = min(ColorBox,1);
else
    ColorBox = ones(PN,1)*ColorFix + [rand(PN,2)*.2-.1 zeros(PN,1)];
    ColorBox = max(ColorBox,0);ColorBox = min(ColorBox,1);
end
  %  hold on
set (gca,'Ydir','reverse')
set (gca,'Xdir','reverse')
set(gca, 'XAxisLocation', 'Top')
set(gca, 'YAxisLocation', 'Right')
PlotI = (randperm(PN))';%
[~, PlotI] = sort(plotZ, 'Descend'); % to plot larger ones first
warning('SORTING!!!')
for px = PlotI'
    loglog(plotX(px),plotY(px), 'o', ...
        'markersize', plotZ(px), ...
        'linewidth', .5, ...
        'MarkerEdgeColor', ColorBox(px,:).^3, ...%[.0 .0 .25], ...
        'MarkerFaceColor', ColorBox(px,:))
hold on
end
box on
grid on
end