% This function makes a pie chart with given labels and colors, up to THR
% minimum frequency. the input "Data" is a verctos of count correponding to
% each thing.
% Pej Oct 2013
function Fig = Pej_Pie_THR(Data, Labels, TypeColor, THR)
if nargin<4
    THR = 3;  % minimum frequency of a Type for being shown in percentage.
end

if nargin<3
    TypeColor = jet(length(Data));
end

tmpCNTr= Data ./ sum(Data) * 100;
tmpFilt = tmpCNTr>THR;

tmpData   = [Data(tmpFilt); sum(Data(~tmpFilt))];
tmpLabels = [Labels(tmpFilt); {'Other'}];

tmpLabels = strrep(tmpLabels, '_', ' ');
Freqs = round(tmpData ./ sum(tmpData) * 1000) / 10;
Freqs(end) = 100 - sum(Freqs(1:end -1)); % just to make sure it sums up to 100.
for i = 1:length(tmpLabels)
    tmpLabels{i} = sprintf('%s\n(%.1f%%)', tmpLabels{i}, Freqs(i));
end

Fig = figure;
h = pie(tmpData, tmpLabels);
set(h(1:2:end), 'EdgeColor', [.991 .991 .991], 'linewidth', 1);
colormap([TypeColor(tmpFilt,:); [.8 .8 .95]]);

end