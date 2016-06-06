function Sim = miRNA_Feat_Filt(miRNA_Data, HealthyIndex, CancerIndex)
global Consts

if Consts.AnalysisMode == 'B'
    Sim = miRNA_Feat_Filt_B(miRNA_Data, HealthyIndex, CancerIndex);
else
    Sim = miRNA_Feat_Filt_C(miRNA_Data, HealthyIndex, CancerIndex);
end
fprintf('All in all %d miRNAs are included in the anlaysis.\n', length(Sim.Data.GenIDs))
end