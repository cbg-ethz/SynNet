% THis file gets a Data matrix containing the couns for different things in
% diffeent samples, and plots them in separate stacked bar plots. 
% Each column in data corresponds to a samples (>SampleLabels), and each row is a single
% Type of thing (>Labels).
% THR is the minimum for maximum frequency in something for being shown


function Fig = Pej_Stacked_Frquencies(Data, Labels, SampleLabels, THR)
if nargin<4
    THR = 3;  % minimum frequency of a Type for being shown in percentage.
end



Frqs= Data ./ repmat(sum(Data,1),size(Data,1),1) * 100;
Filt = max(Frqs,[],2)>THR;
FrqsF = Frqs(Filt,:);
LabelsF = Labels(Filt);
FrqsF(end+1,:) = 100 - sum(FrqsF,1);
LabelsF(end+1) = {'Other'};
LabelsF = strrep(LabelsF, '_', ' ');
Fig = figure;
CBox = jet(size(FrqsF,1)); CBox = CBox(randperm(size(FrqsF,1)),:);
colormap(CBox);
h = bar(-FrqsF', 'BarLayout', 'Stacked', 'Barwidth', .9); % I inverse them as a quick hack cuz matlab automakes thislagend otherway around
%set(h, 'EdgeColor', [.9 .9 .9]);
legend(LabelsF, 'Location','EastOutside');legend('boxoff'); 
set(gca, 'XTick', 1:length(SampleLabels));
set(gca, 'XTickLabel', SampleLabels);
ylim([-102.5 2.5])
set(gca, 'YTick', -100:20:0);
set(gca, 'YTickLabel', 0:20:100);
xlim([0 length(SampleLabels)]+.5)
box off
set(gcf, 'position', [680         868        1000         250]);

end
