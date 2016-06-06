function Pej_Test_DEG_Enrichments(DEGoutputfolder,DB_Path)
Qthr = 0.01;
Shuffle = false; % If you put this on true, it shuffles the DEG qvalues, so should technically give flat pvalues all the time.


if nargin < 2
    DB_Path = '/Users/pejmanm/Desktop/LocalTMP/PEJ_Resources/GeneSets/';
end

Resultfiles = dir(DEGoutputfolder);
% Pej_GetFiles([DEGoutputfolder '/*.txt']);
% Resultfiles = regexp(Resultfiles, '[\f\n\r]', 'split');
for f = 1 : length(Resultfiles)
    if ~isempty(Resultfiles)
       [~, Fname, Fext] =fileparts(Resultfiles(f).name);
        if ~strcmpi(Fext, '.txt')
            continue
        end
        disp(['Analysing ' Resultfiles(f).name])
        Test_DEGs_Enrichments_of(DEGoutputfolder, Fname,DB_Path, Shuffle, Qthr);
    end
end
end


function Test_DEGs_Enrichments_of(Fldr, Fname,DB_Path, Shuffle, Qthr)
% if iscell(DESeq_Output); DESeq_Output=DESeq_Output{1};end

% [Fldr, Fname, Fext] =fileparts(DESeq_Output);
OutFldr = [ Fldr '/DEG_EnrichmentTests'];
mkdir(OutFldr)
Clustering_Output = [OutFldr '/Clusters_' Fname];
DES_Res = Pej_Read_Table([Fldr '/' Fname '.txt']);
DES_Res = Pej_Struct_RowDel(DES_Res,isnan(DES_Res.padj));


%% Clustering
Clustered_Genes.IDs = strrep(DES_Res.UnLabeled_C1, '"', ''); % the background for the enrichments
if Shuffle
Qvals = DES_Res.padj(    randperm(length(DES_Res.padj))); % Shuffle
else
Qvals = DES_Res.padj; % Qvas of test
end
CF    = DES_Res.log2FoldChange;

Clustered_Genes.ClusterLabels2{1} = ['Upregulated'];
IDXKM(:,1)= Qvals<Qthr & CF>0;

Clustered_Genes.ClusterLabels2{2} = ['Dnregulated'];
IDXKM(:,2)= Qvals<Qthr & CF<0;

Clustered_Genes.Cluster_IDs     = IDXKM;
Clustered_Genes.EffectSize      = CF;

save(Clustering_Output, 'Clustered_Genes');

Pej_Test4Enrichment(Clustering_Output, DB_Path);
clear Clustered_Genes IDXKM
end