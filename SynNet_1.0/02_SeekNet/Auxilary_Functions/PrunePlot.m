function PrunePlot(T)
% T = Pej_Struct_RowSelect(T, T.Samples<=100);
figure;
hold on
subplot(1,2,1)
plot([50 100], [50 100], 'g')
hold on
plot(T.AUC_L, T.AUC_Lpru, '.')

subplot(1,2,2)
plot([50 100], [50 100], 'g')
hold on
plot(T.AUC_T, T.AUC_Tpru, '.')

Recall    = [T.Syn_and_Best, T.Syn_and_Prun]./ [T.GenesinSynthCirc  T.GenesinSynthCirc ];
Precision = [T.Syn_and_Best, T.Syn_and_Prun]./ [T.GenesinLearntCirc T.GenesinPrunedCirc];

figure
hold on
plot(Recall(:,2), Precision(:,2), 'o')
plot(Recall(:,1), Precision(:,1), '.r')
end