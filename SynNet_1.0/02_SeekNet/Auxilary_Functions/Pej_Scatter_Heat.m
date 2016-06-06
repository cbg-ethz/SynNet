function Pej_Scatter_Heat(plotX,plotY,plotZ, ColorMap)
PN = size(plotX,1);
if nargin < 4
    ColorMap = 'jet';
end
Nbins = 40;

Cbox = eval([ColorMap '(ceil(Nbins*1.2));']);
Cbox = Cbox(round(Nbins/10): round(Nbins/10)+Nbins-1,:); % remove two extreme of the color.

Zqs = quantile(plotZ, linspace(0,1, Nbins+1));
Zqs(1) = [];
CI = nan(size(plotZ));
for i = Nbins:-1:1
    CI(plotZ<=Zqs(i)) = i;
end
hold on
set (gca,'Ydir','reverse')
set (gca,'Xdir','reverse')
set(gca, 'XAxisLocation', 'Top')
set(gca, 'YAxisLocation', 'Right')

PlotI = randperm(PN); % plot points in an arbitrary order, in order to make sure there won't be any artificial collor patterns dues to the plot orderring.
for px = PlotI
    plot(plotX(px),plotY(px), 'o', ...
        'markersize', 4, ...
        'linewidth', 1, ...
        'MarkerEdgeColor', Cbox(CI(px),:), ...
        'MarkerFaceColor', 'none')%Cbox(CI(px),:))
end
box on



end