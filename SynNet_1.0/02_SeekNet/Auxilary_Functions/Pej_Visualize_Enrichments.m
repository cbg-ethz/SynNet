function Pej_Visualize_Enrichments(OutputFolder, EnrichmentOutputPathPattern, varargin)
Q_Thr = .01;
MinQ_log10 = 1E-15; % The smallest Qvalue to plot, smaller ones saturate here.

mkdir(OutputFolder);

inFiles = Pej_GetFiles(EnrichmentOutputPathPattern);
NinF = length(inFiles);
fprintf('%d input files found.\n', NinF);

if ~isempty(varargin)
    for j = 2:length(varargin)+1
        inFiles(j,:) = Pej_GetFiles(varargin{j-1});
    end
end

for i = NinF:-1:1
    [~,Fname{i},~] = fileparts(inFiles{1,i});
    D(i)    = load(inFiles{1,i});
    for j = 2:size(inFiles,1)
        DH      = load(inFiles{j,i});
        %         for hset = 1:length(DH.Input.GeneSets)
        %             DH.Input.GeneSets(hset).Genes = [DH.Input.GeneSets(hset).Genes -1]; % I add the invalid genes ID "-1" to all HIV related pathways to make sure they won't become a subset of anything else.
        %         end
        D(i).Input.GeneSets           = [D(i).Input.GeneSets             DH.Input.GeneSets];
        D(i).Results.meancoeffs       = [D(i).Results.meancoeffs  ; DH.Results.meancoeffs];
        D(i).Results.Q_val            = [D(i).Results.Q_val            ; DH.Results.Q_val];
    end
    Sigs(:,:,i) = D(i).Results.Q_val <= Q_Thr;
end
%% Remove Duplicate sets (This can particularely happen if more than one DB is used)
IsuniqGS = Get_Uniques(D);

All_Sigs = find(IsuniqGS & any(any(Sigs,3),2)); clear Sigs

GS = D(1).Input.GeneSets;
for i = length(All_Sigs):-1:1
    SGS.Ls(i) = length(GS{All_Sigs(i)}.Genes);
    SGS.Names{i} = GS{All_Sigs(i)}.Name{1};
    SGS.Genelist{i} = GS{All_Sigs(i)}.Genes;
end
SGS.Names{end+1} = 'Miscellaneous';

%% Find Main Categories

SGS.SuperSet = WhosYourDaddy(SGS.Genelist);
[SGS.X_cords,  BorderCords] = Distribute(SGS.SuperSet);
Dads = unique(SGS.SuperSet);
%% Prepare colors
ColorBox(Dads,:)  = lines(length(Dads));
% Manually change some of the colors!
ColorBox(Dads(  end  ),:)  = [.7 .7 .7]; % The 'Miscellaneous' group

for  j = 1:length(ColorBox)-1
    ColorBox(j,:) = ColorBox(SGS.SuperSet(j),:)+rand(1,3)*0.1;%.*(1+randn(1,3)*0.1);
end
ColorBox = max(ColorBox,0);ColorBox = min(ColorBox,1);
%% Plot!
for i = 1:NinF
    tmpD = D(i);
    %   tmpC = C(i);
    iSGS.Qs = tmpD.Results.Q_val(All_Sigs,:);
    iSGS.AVGeffect = tmpD.Results.meancoeffs(All_Sigs,:);
    
    %% Push data in bounds for good plotting
    iSGS.Qs(~isnan(iSGS.Qs)) = max(iSGS.Qs(~isnan(iSGS.Qs)),MinQ_log10);
    
    iSGS.AVGeffect(iSGS.AVGeffect> 3) =  3;
    iSGS.AVGeffect(iSGS.AVGeffect<-3) = -3;
    
    %% plotting
    Fig = figure;
    set(gcf, 'position', [1 1 1000 250])
    hold on
    plot([0.025 max(BorderCords)-.05], [0 0], '-', 'linewidth', 1, 'color', [0 0 .2])
    %title(Stimuli{i})
    box off
    grid off
    for j = 1:length(SGS.X_cords)
        if any(~isnan(iSGS.Qs(j,:)))
            %% This is ad-hoc edited for DEG results, does not go with more than two clusters
            %[tmpQ tmpI] = min(iSGS.Qs(j,:));
            tmpQ = iSGS.Qs(j,1);tmpI = 1;
            
            if tmpQ<=Q_Thr
                plot(SGS.X_cords(j),iSGS.AVGeffect(j, tmpI), 'o', ...
                    'markersize', 2+(-log10(tmpQ)), ...
                    'linewidth', 1, ...
                    'MarkerEdgeColor', ColorBox(SGS.SuperSet(j),:), ...
                    'MarkerFaceColor', ColorBox(j,:))
            end
            
            tmpQ = iSGS.Qs(j,2);tmpI = 2;
            
            if tmpQ<=Q_Thr
                plot(SGS.X_cords(j),iSGS.AVGeffect(j, tmpI), 'o', ...
                    'markersize', 2+(-log10(tmpQ)), ...
                    'linewidth', 1, ...
                    'MarkerEdgeColor', ColorBox(SGS.SuperSet(j),:), ...
                    'MarkerFaceColor', ColorBox(j,:))
            end
            
            
        end
    end
    set(gca, 'XTick', []);
    xlim([-.025 max(BorderCords)])
    ylim([-3.25 3.25])
    set(gca, 'YTick', [-3:3]);
    set(gca, 'XTick', []);
    %ylabel('Fold change(log_2)')
    for j = 1:length(Dads)-1
        plot(BorderCords(j)*[1 1], ylim *.9 , '--','color', [0 0 .2],'linewidth', 1)
    end
    
    % Add my own axes
    % Leftside
    plot(zeros(1,7)-0.025, -3:3, 'k+','linewidth', 1)
    plot([0 0]-0.025, [-3.25 +3.25], 'k','linewidth', 1.2)
    text(zeros(1,7)-0.025 -.01, -3:3, mat2cell(-3:3), 'HorizontalAlignment', 'right', 'FontSize', 14)
    % Right Side
    plot(ones(1,7)*BorderCords(end), -3:3, 'k+','linewidth', 1)
    plot([1 1]    *BorderCords(end), [-3.25 +3.25], 'k','linewidth', 1.2)
    text(ones(1,7)*BorderCords(end)+0.05, -3:3, mat2cell(-3:3), 'HorizontalAlignment', 'right', 'FontSize', 14)
    
    
    
    axis off
    Pej_SavePlot(Fig, [OutputFolder '/' Fname{i}]);
end

%% Make Legened plot
Fig = figure; hold on
XD = [-0.025 BorderCords];
for j = 1:length(Dads)
    jXD = mean(XD(j:j+1));
    plot(jXD, 0, 'o', ...
        'markersize', 15, ...
        'linewidth', 1, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', ColorBox(Dads(j),:));
    
    DN = strrep([SGS.Names{Dads(j)} '  '], '_', ' ');
    text(jXD, 0, DN, 'rotation', 45, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14);
    %plot(BorderCords(j)*[1 1], ylim, '--', 'color', [.8 .8 .8])
    xlim([-.025 max(BorderCords)])
    
    set(gcf, 'position', [1 1 1000 250])
    axis off
end
Pej_SavePlot(Fig, [OutputFolder '/Legend']);
%% Make SizeLegend plot
Fig = figure; hold on
for p = 5:5:15
    plot(0, p/5, 'o', ...
        'markersize', 2+p, ...
        'linewidth' , 1, ...
        'MarkerEdgeColor', [0.3 0.3 0.3], ...
        'MarkerFaceColor', [0.2 0.2 0.3]);
    text(.05, p/5, ['10^{', int2str(-p), '}'], 'HorizontalAlignment', 'left' );
    
    set(gcf, 'position', [1 1 1000 250])
    xlim([-.025 max(BorderCords)])
    ylim([-3.25 3.25])
    axis off
end
Pej_SavePlot(Fig, [OutputFolder '/SizeLegend']);

end
function [X_cords, BorderCords] = Distribute(SuperSets)
X_cords = nan(size(SuperSets));
Dads = unique(SuperSets);
M = length(SuperSets);
N = length(Dads);
X0 = 0;dX = .025;
for d = 1:length(Dads)
    Kids = SuperSets==Dads(d);
    KL = sum(Kids) / M;
    X_cords(Kids) = X0 + dX + KL * rand(sum(Kids),1);
    BorderCords(d)   = X0 + dX + KL + dX;
    X0 = X0 + dX + KL + dX;
end
end
function [Daddy_List] = WhosYourDaddy(Genelists, Ndaddies)
THR = .95; % Percentage for Chilhood threshold
if nargin < 2
    Ndaddies = 9;
end

N = length(Genelists);
Intsct = zeros(N);
for i = 1:N-1
    Si = Genelists{i};
    for j = i+1:N
        Sj = Genelists{j};
        intL = length(intersect(Si,Sj));
        Intsct(i,j) = intL / length(Si);
        Intsct(j,i) = intL / length(Sj);
    end
end

Maximals = max(Intsct,[],2)<THR; % These guys are less than 95% redundant with somemthing else.
KidSum   = sum(Intsct,1)'; % how much each set inculdes the others.
[~, Daddies] = sort(KidSum .* Maximals, 'Descend');

DaddiesFilt = false(N,1);
DaddiesFilt(Daddies(1:Ndaddies)) = true;
Intsct(:,~DaddiesFilt) = 0;
[MaxInt, Daddy_List] = max(Intsct+rand(size(Intsct))*1E-3, [], 2);
Daddy_List(MaxInt<=THR) = N+1;
Daddy_List(DaddiesFilt) = find(DaddiesFilt);
end


function IsUniq = Get_Uniques(D)
% look into the first file set, and find duplicates!
nD = length(D);
for i = length(D(1).Input.GeneSets):-1:1
    lGS(i,1) = length(D(1).Input.GeneSets{1,i}.Genes_PejIDs);
    pGS(i,1) = 1;
    for j = nD:-1:1
        pGS(i,1) = min([pGS(i), D(j).Results.Q_val(i,:)]);
    end
end

[sL, sLI] = sortrows([lGS, pGS], [1 1]); % here I sort also by the best p-value so that I discard the duplicates with worse pvalues always
IsUniq = true(size(pGS));
for i  = 1:length(lGS)-1
    if sL(i)~=sL(i+1)
        IsUniq(sLI(i+1)) = true;
    else
        % check if they have equal genes
        if sL(i) == length(intersect(D(1).Input.GeneSets{1,sLI(i)}.Genes_PejIDs, D(1).Input.GeneSets{1,sLI(i+1)}.Genes_PejIDs))
            % they are equal
            IsUniq(sLI(i+1)) = false;
        end
        
    end 
end


end