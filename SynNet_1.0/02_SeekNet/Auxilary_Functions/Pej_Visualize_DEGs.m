function Pej_Visualize_DEGs(DEGoutputfolder, FixPlotLims)
Qthr = 0.01;
TextLabel_Thr = 0; %  Genes with q-value smaller than this will be plotted with their names printed next to them.

if nargin<2
    FixPlotLims=false; % make plots independently
end
%--------
Resultfiles = dir(DEGoutputfolder);
MinY = 1;
MaxX = 0;
if FixPlotLims
    for f = 1 : length(Resultfiles)
        if ~isempty(Resultfiles)
            [~, Fname, Fext] =fileparts(Resultfiles(f).name);
            if ~strcmpi(Fext, '.txt')
                continue
            end
            tmpDEseq = Pej_Read_Table(Resultfiles(f).name, [], false);
            MinY = min([MinY; tmpDEseq.padj(tmpDEseq.padj>0)]);
            MaxX = max([MaxX; abs(tmpDEseq.log2FoldChange)]);
        end
    end
end

for f = 1 : length(Resultfiles)
    if ~isempty(Resultfiles)
        [~, Fname, Fext] =fileparts(Resultfiles(f).name);
        if ~strcmpi(Fext, '.txt')
            continue
        end
        disp(['Analysing ' Resultfiles(f).name])
        Get_DEG_stats(Resultfiles(f).name, Qthr, f==1);
        PlotForOne(DEGoutputfolder, Fname, Qthr, TextLabel_Thr, MinY, MaxX);
    end
end
end


function PlotForOne(Fldr, Fname, Qthr, TextLabel_Thr, MinY, MaxX)
% set(0, 'defaultTextInterpreter', 'latex');
OutFldr = [ Fldr '/DEG_Plots'];
mkdir(OutFldr)
DES_Res = Pej_Read_Table([Fldr '/' Fname '.txt'], [], false);
DES_Res = Pej_Struct_RowDel(DES_Res, isnan(DES_Res.log2FoldChange));

% Polish the qvalues a bit
DES_Res.padj(isnan(DES_Res.padj))=1; % give 1 to those that were excluded from test
if any(DES_Res.padj==0)
    newZero = min(DES_Res.padj(DES_Res.padj>0))/10;
    DES_Res.padj(DES_Res.padj==0) = newZero; % put a small value for those that hit machine's epsilon
else
    newZero=0;
end

% plot it!
Fig = figure;
plotX = DES_Res.log2FoldChange;
plotY = max(-inf, DES_Res.padj);
plotZ = log(1+DES_Res.baseMean);
plotL = strrep(DES_Res.UnLabeled_C1, '"', '');
PN = size(plotX,1);
Plot_H1 = min(1, sum(DES_Res.padj>Qthr)/1000); % Probability of of ploting a non-significant point. I discard some to avoid making too heavy plots

ColorBox = ones(PN,1)*[.15 .15 .7] + [rand(PN,2)*.3-.15 zeros(PN,1)];
ColorBox = max(ColorBox,0);ColorBox = min(ColorBox,1);
semilogy(xlim, [1 1]* Qthr, 'r')
hold on
if MinY<1 && MaxX>0
    xlim(ceil(MaxX)*[-1 +1]);
    ylim([MinY*1E-5 10^(-log10(MinY)/40) ])
end
[~, PlotI] = sort(plotZ, 'Descend'); % to plot larger ones first
for px = PlotI'
    if plotY(px)>Qthr && rand < Plot_H1
        % skip it!
        continue
    end
    semilogy(plotX(px),plotY(px), 'o', ...
        'markersize', plotZ(px), ...
        'linewidth', min(1,plotZ(px)*.1), ...
        'MarkerEdgeColor', [.0 .0 .25], ...
        'MarkerFaceColor', ColorBox(px,:))
    if plotY(px)<TextLabel_Thr
        text(plotX(px),plotY(px),plotL{px},'fontsize', 8,'verticalalignment', 'top')
    end
end
xlim([-1 1]* (max(abs(xlim))));
% ylim([1E-13 2]);
% set(gca, 'YTick', 10.^[-12:2:0])
semilogy(xlim, [1 1]* Qthr, 'r')
box on
% text(min(xlim)+.1, 0.06, sprintf(' %.0f\% FDR',Qthr*100), 'color', [.3 0 0], 'verticalalignment', 'bottom', 'horizontalalignment', 'left','fontsize', 9, 'fontweight', 'bold');
xlabel('Fold Change (log_2)')
ylabel('Adjusted p-value')
set(gcf, 'position', [1 1 330 240]);

%% Make SizeLegend plot
cnt = .5; MY = 0;%-min(log10(ylim));
dy  = range(log10(ylim)).*.1;
dx  = range(xlim).*.1;
for SIZE = logspace(1,4, 4)
    p = log(1+SIZE);
    semilogy(min(xlim)+dx*.3, 10^-(MY+cnt*dy), 'o', ...
        'markersize', p, ...
        'linewidth' , min(1,p*.1), ...
        'MarkerEdgeColor', [.0 .0 .2], ...
        'MarkerFaceColor', [0.7 0.7 0.7]);
    text(min(xlim)+dx*.7, 10^-(MY+cnt*dy), ['10^{', int2str(log10(SIZE)), '}'], 'HorizontalAlignment', 'left','fontsize', 9);%, 'fontweight', 'bold' );
    
    cnt = cnt+.5;
end
ylim([min(ylim) 10^(dy/4)])

XL = xlim;
YL = ylim;
fill(XL([1 2 2 1 1]), [1 1 YL([2, 2]) 1], 'w')

%% Print the number of significant and non-sig genes
text(max(xlim)-dx*.2, 10^-(MY+0.5*dy), sprintf('H0:%d', sum(DES_Res.padj>Qthr)), 'HorizontalAlignment', 'right','fontsize', 9);%, 'fontweight', 'bold' );
text(max(xlim)-dx*.2, 10^-(MY+1.0*dy), sprintf('H1:%d', sum(DES_Res.padj<=Qthr)), 'HorizontalAlignment', 'right','fontsize', 9);%, 'fontweight', 'bold' );


PejSavePlot(Fig, [OutFldr '/' Fname])
% set(0, 'defaultTextInterpreter', 'tex');
end

function PejSavePlot(Fig, OutFig, Boxoff)
if nargin <3; Boxoff = false;end

saveas(Fig, OutFig)
if Boxoff;
    h = get(Fig,'Children');
    for i = 1:length(h)
        try
            subplot(h(i));
            legend boxoff;
        catch ignore
            
        end
    end
end

set(Fig, 'PaperPositionMode', 'auto');
print(Fig, '-depsc2', OutFig);
close(Fig)
disp(['Figure saved: ' OutFig] )
end

function Get_DEG_stats(File, QThr, Append)
% QThr = .01;
if Append
    Fout = fopen('DEGstats.log', 'a');
else
    Fout = fopen('DEGstats.log', 'w');
    fprintf(Fout, 'ComparisonID\tDEGs(%d%% FDR)\tDEGs&2foldchange\tTotalTested\n', round(QThr*100));
end
F = Pej_Read_Table(File, [], false);
fprintf(Fout, '%s\t%d\t%d\t%d\n', File,  sum(F.padj<=QThr), ...
    sum(F.padj<=QThr & abs(F.log2FoldChange)>1), sum(~isnan(F.padj)));
fclose(Fout);

end