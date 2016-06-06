% NOTE: This code is not yet matured. 
% Tuning the layout and size of the window is still independant of the data
% dimensions and it's manual.



% if you want the arrows not to collide sort the data in descending order
% based on the rowwise means:
% [~, I] = sort(mean(Data,2), 'descend');
% Pej_Heatmap(Data(I,:));
 

% Pej 2014, Sept.
function [Fig, h] = Pej_Heatmap(Data, RowNames, ColumnNames, Scale)
N = size(Data,1);
P = size(Data,2);
CRange = [min(Data(:)) max(Data(:))];%[min(Top_miRNAs.RelativeFreq_log10(:)), max(Top_miRNAs.RelativeFreq_log10(:))];
BM = mean(Data,2);

if nargin<2 || isempty(RowNames)
    RowNames = cell(N,1);
end

if nargin<3 || isempty(ColumnNames)
    ColumnNames = cell(P,1);
end

if nargin<4
    Scale = 'absolute';
end
switch Scale
    case 'log2'
       ScalePrefix = '2^';
    case 'log2'
       ScalePrefix = '10^';
    otherwise
        ScalePrefix = '';
end
%% initialize the figure layout
Fig = figure;
set(gcf, 'position', [500 500 400 650])
h(1) = subplot('Position',[.1 .1 .15 .8]);
h(2) = subplot('Position',[.25 .1 .4 .8]);
h(3) = subplot('Position',[.65 .1 .07 .8]);
h(4) = subplot('Position',[.72 .1 .06 .8]);
h(5) = subplot('Position',[.8 .1 .09 .8]);

%% Make the heatmap
subplot(h(2));
Cb = bone(255).^1.5;
Cb(end:-1:1,:)= Cb;
colormap(Cb);

imagesc(Data)

set(gca, 'YDir', 'reverse');
set(gca, 'YTick',[]);
set(gca, 'YTickLabel', []);
set(gca, 'XTick', []);
YL = ylim;
YLc = [1 size(Cb,1)];
tmph = colorbar;
set(tmph, 'Ylim', CRange);
MainYticks = get(tmph, 'YTick');

colorbar('off')


%% Make names
subplot(h(1));
text(zeros(N,1), .25+[N:-1:1], RowNames, 'fontname', 'calibri', 'horizontalalignment', 'left', 'verticalalignment', 'top', 'fontsize', 9)
ylim(YL);
xlim([0 1])
axis off
%% Make arrows
subplot(h(3));
axis off
caxis(CRange);
% Regress the Basemeans to the color scale
BMc = Expr2ColorIndex(BM, CRange, YLc);
BMi = round(BMc);
Pos = 1:N;
Posc = Expr2ColorIndex(N-Pos+1, YL, YLc);
hold on
for i = 1:length(BM)
    plot([0 1], [Posc(i) BMc(i)], 'color', Cb(BMi(i),:), 'linewidth', 1);
end
ylim(YLc);
xlim([0 1])
axis off
%% make colorbar
subplot(h(4));
imagesc(linspace(CRange(1), CRange(2), 255)')
set(gca, 'YAxisLocation', 'right');
% Yticks = log10(1./[1:9 10:10:90 100:100:1000]); Yticks = Yticks(end:-1:1);
Yticks = MainYticks;
% Regress the Ticks to the color scale
Yticks = Expr2ColorIndex(Yticks, CRange, YLc);

set(gca, 'YTick', Yticks);
set(gca, 'XTick', []);
set(gca, 'YTickLabel', []);
set(gca, 'YDir', 'normal');

ylim(YLc);
colormap(Cb);
caxis(CRange);
hold on
%% Make ColorBar legend
subplot(h(5));
ylim(YLc);
set(gca, 'YDir', 'normal');

axis off
for i = 1:length(MainYticks)
    text(0, Expr2ColorIndex(MainYticks(i), CRange, YLc), [ScalePrefix '{' num2str(MainYticks(i)) '}'],'fontname', 'calibri', 'horizontalalignment', 'left', 'verticalalignment', 'middle')
end
xlim([0 1])
end


function ColorI = Expr2ColorIndex(Exp, Erange, Crange)
% This regresses between color Index and expressions
ColorI = (Exp - min(Erange))/range(Erange)*range(Crange)+min(Crange);


end