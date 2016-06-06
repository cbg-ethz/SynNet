function Write_Function_InOuts_B(Function, Data, Results, FigurePath)

Max_Range = 2;
WindowWidth = 650;
Function = squeeze(Function);
Function = (sortrows(Function', 1))'; % make sure if there's something in a row it's in the begining
Function = sortrows(Function, 1); % Sort such that NOT branches are first
if numel(Data.Threshold)==1; Data.Threshold=Data.Threshold*ones(size(Data.Values));end

if ~any(Function(:))
    % Its an empty function
    Fig = figure;
    plot(1,1);
    text(1, 1,'Empty Circuit!');
    Pej_SavePlot(Fig, FigurePath)
    return
end


A = Data.Annots;
BranchSep = [];Genes =[];
for i = 1: size(Function,1)
    Fr = Function(i,:)~=0;
    if any(Fr)
        Genes = [Genes abs(Function(i,Fr))];
        BranchSep(end+1) = length(Genes)+.5;
    end
end
GenesL= strrep(lower(Data.GenIDs(Genes))        , 'hsa-', '');


Ns = length(A);
Ng = length(Genes);
DotSZ = min(5, WindowWidth*.8/Ns/2);
% D = Data.Values(Genes,:) - log2(Data.Threshold(Genes,:));
% Os = nan(Ns,1);
% if Ng>1
%     % Order samples by hierarchical clustering
%     if sum(A)>1
%         Os_p = clusterdata(D(:, A)','cutoff', eps, 'distance', 'cosine');
%     else
%         Os_p = 1;
%     end
%     if sum(~A)>1
%         Os_n = clusterdata(D(:,~A)','cutoff', eps, 'distance', 'cosine');
%     else
%         Os_n=1;
%     end
%     Os( A) = Os_p + sum(~A);
%     Os(~A) = Os_n;
% else
%     % Order samples by expression
%     Os( A) = D( A)/min(D( A))+eps;
%     Os(~A) = D(~A)/max(D(~A))-eps;
%     
% end
% 
% [~, SampleOrder] = sort(Os, 'Ascend');
TordM = Results.Output_B_Margins;
TordM(~Results.Output_B) = -TordM(~Results.Output_B);
[~, SampleOrder] = sortrows([double(A); TordM]', [1 2]);
[~, SampleOrder] = sortrows([TordM]', [ 1]);
 
D2plot = Data.Values (Genes,SampleOrder);
B2plot = Data.BValues(Genes,SampleOrder);

D2plot = log2(D2plot);
D2plot = D2plot - log2(Data.Threshold(Genes,SampleOrder));
D2plot = D2plot/Max_Range;
CRange = [-1 1];
XLs = [.5, Ns+.5];
YLs = [.5, Ng+.5];

WindowHeight = 100+(Ng+2.75)*20;
TextPart = (  50/WindowHeight) ;
GraphPart= (1-50/WindowHeight) ;
PlotSize = GraphPart / (Ng + 1 + .5 + 1 + .25);

Fig = figure;
set(gcf, 'position', [500 500 WindowWidth  WindowHeight])

%% Function Inputs
subplot( 'position', [.2 PlotSize*2.75 .8 PlotSize*Ng])
colormap(jet)
imagesc(D2plot)
hold on
caxis(CRange);
for i = 1:size(B2plot,1)
    tmpF = find(B2plot(i,:));
    if~isempty(tmpF)
        plot(tmpF, i, 's', 'markersize', DotSZ, 'markerfacecolor', [0.9 0 .1], 'color', [.01 0 0.1 ], 'linewidth', min(.5, DotSZ/2));
    end
end
set(gca, 'XTick', []);
set(gca, 'YTick', []);
if Ng>1
    plot(xlim, [1;1] * (1.5:(Ng-.5)), '-' , 'linewidth', .5, 'color', [.8 .8 .8 ])
end
BranchSep(BranchSep>Ng)=[];
if ~isempty(BranchSep)
    plot(xlim, [1;1] * BranchSep , '-', 'linewidth', .5, 'color', [0  0  .01])
end
xlim(XLs);
ylim(YLs);
Ydir = get(gca, 'Ydir'); % This is crucial! Because different matlab versions do this differently, if you don't explicitely fix the direction of your axis you might end up having a reverse label for your genes.
Xdir = get(gca, 'Xdir'); % This is crucial! Because different matlab versions do this differently, if you don't explicitely fix the direction of your axis you might end up having a reverse label for your genes.


%% Sample Labels
subplot( 'position', [.2 GraphPart .8 TextPart ])
set(gca, 'Xdir', Xdir);
for i = 1:Ns
    text(i, .1,Data.SampleID{SampleOrder(i)}, 'horizontalalignment', 'left', 'rotation', 90, 'fontsize', DotSZ*2)
end
axis off
xlim(XLs);
ylim([0 1]);
%% Input Labels
subplot( 'position', [0 PlotSize*2.75 .2 PlotSize*Ng])
set(gca, 'Ydir', Ydir);
for i = 1:Ng
    text(.95,i, GenesL{i}, 'horizontalalignment', 'right')
end
axis off
xlim([-1 1])
ylim(YLs);

%% Function output
subplot( 'position', [.2 PlotSize*1.25 .8 PlotSize ])
CO = Results.Output_B_Margins*log2(10); % in log2
CO(~Results.Output_B)=-CO(~Results.Output_B);
COm = 0;% (mean(CO(A)) + mean(CO(~A)))/2;
CO  = CO - COm;
CO  = CO / Max_Range;


imagesc(CO(SampleOrder));
caxis(CRange);
hold on
if any(Results.Output_B)
    plot(find(Results.Output_B(SampleOrder)), 1, 's', 'markersize', DotSZ, 'markerfacecolor', [0.9 0 .1], 'color', [.01 0 0.1 ], 'linewidth', min(.5, DotSZ/2));
end

set(gca, 'XTick', []);
set(gca, 'YTick', []);
xlim(XLs);
ylim([.5 1.5]);


%% Output Label
subplot( 'position', [0 PlotSize*1.25 .2 PlotSize])
text(.95,1, 'Output', 'horizontalalignment', 'right')
axis off
xlim([-1 1])
ylim([.5 1.5]);

%% True Annotation
subplot(  'position', [.2 0 .8 PlotSize ])
imagesc((A(SampleOrder)-.5)*2);
% caxis(CRange);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
xlim(XLs);
ylim([.5 1.5]);

%% True Label
subplot( 'position', [0 0 .2 PlotSize])
text(.95,1, 'Truth', 'horizontalalignment', 'right')
axis off
xlim([-1 1])
ylim([.5 1.5]);

Pej_SavePlot(Fig, FigurePath)

%% Color Legend
Fig = figure;
set(gcf, 'position', [500 500  WindowWidth (Ng+2.5)*25]);
colormap(jet)
colorbar
caxis(CRange* Max_Range);
axis off
Pej_SavePlot(Fig, [FigurePath '-Clegend'])

%% Make a report file


end