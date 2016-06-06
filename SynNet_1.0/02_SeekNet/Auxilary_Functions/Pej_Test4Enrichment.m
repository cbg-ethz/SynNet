% Before using this you should once have imported the GMT files from the
% DBPath using: "Pej_ImportGMT".
function OutPutPrefix = Pej_Test4Enrichment(Clustering_Output, DB_Path, XrefPath)
Thr = 0.05; % default FDR threshold

if nargin < 2
    DB_Path = '/Users/pejmanm/Desktop/LocalTMP/PEJ_Resources/GeneSets/';
end

if nargin < 3
    XrefPath = '/Users/pejmanm/Desktop/LocalTMP/PEJ_Resources/Gene-Xref-25Apr2015.txt';
end

listing = dir(DB_Path);
ValidDB = false(length(listing),1);
for dbFile = 1:length(listing)
    [Fpath FName Fext] = fileparts(listing(dbFile).name);
    ReferenceFilesLabels{dbFile} = FName;
    ReferenceFiles{dbFile} = [DB_Path '/' FName Fext];
    
    if isempty(FName) || FName(1) == '.' || ~strcmpi(Fext, '.mat'); continue; end
    
    tmpdb = load(ReferenceFiles{dbFile});
    if sum(isfield(tmpdb, {'GeneSets' 'GeneNames'}))<2 ; continue; end
    ValidDB(dbFile) = true;
end
ReferenceFilesLabels(~ValidDB)= [];
ReferenceFiles(~ValidDB)= [];

ClPath = Clustering_Output(1:find(Clustering_Output=='/',1,'last'));
ClPfix = Clustering_Output(find(Clustering_Output=='_',1,'last'):end);

OutFldr = [ClPath 'Test4Enrichment' ClPfix];
mkdir(OutFldr);
OverallReport = [Clustering_Output '_Test4Enrichment_Report.txt'];
F = fopen(OverallReport, 'w');
fprintf(F, 'DataBase\tClusterName\tSetName\tOver/UnderRep.\tlog10(Pval)\tlog10(Qval)\tCountinBackground\tExpectedCount\tCountinCluster\tAvgEffectSize\tSourceUrl\tClusterID\ttheGenes\n');

load(Clustering_Output);

for SourceDB = 1:length(ReferenceFilesLabels)
    clear P1 P2 P3 P4
    load(ReferenceFiles{SourceDB});    
    Clustered_Genes = Get_PejIDs(Clustered_Genes, GeneNames,XrefPath);
    P2 = length(Clustered_Genes.Cluster_IDs); % P2 = 2nd parameter in hygecdf(X,M,K,N)
        disp(['Enrichments are caclulated against ' int2str(P2) ' background genes'])

    if isvector(Clustered_Genes.Cluster_IDs)
        % IDX is a vector containing unique cluster IDS, as given by Matlab kmeans function.
        for k = Clustered_Genes.Kchoice:-1:1
            Clusters(:,k) = Clustered_Genes.Cluster_IDs==k;
        end
    else
        % IDX is a membership matrix, whith each column corresponding to one cluster.
        Clusters = Clustered_Genes.Cluster_IDs;
    end
       
    fprintf('%s was loaded for analysis\n', ReferenceFilesLabels{SourceDB});
    minSetSize = 2; % Smallest set/reference set to be considered.
    
    N = length(GeneSets) ; % Total number of background Sets    
    for s = N:-1:1
        P3(s,1) = length(intersect(GeneSets{s}.Genes_PejIDs,Clustered_Genes.Pej_ID));% P3 = 3rd parameter in hygecdf(X,M,K,N)
    end
    
    K = size(Clusters,2); % total number of tested clusters
    NK = N * K; NK10th = ceil(NK/10);
    
    P_val = nan(N, K); SetCounts = nan(N, K);Ex_P1 = nan(N, K);
    clear GenesInClust
    GenesInClust{N,K} = [];
    for k = K:-1:1        
        P4 = sum(Clusters(:,k)); % P4 = 4th parameter in hygecdf(X,M,K,N)
        if P4 <= minSetSize
            % screw it!
            sNk = 0+N*k;
            if mod(sNk,NK10th)==0 && sNk/NK10th < 11
                fprintf('%d%% ',  (sNk/NK10th) * 10);
            end
            continue
        end
        
        
        for s = N:-1:1
            if P3(s)<minSetSize
                continue;
            end
            GenesInClust{s,k} = intersect(Clustered_Genes.Pej_ID(Clusters(:,k)), GeneSets{s}.Genes);
            P1 = length(GenesInClust{s,k});% P1 = 1st parameter in hygecdf(X,M,K,N)
            if P1 <= minSetSize
                P_val(s, k) = nan;
            else
                P_val(s, k) = hygecdf(P1-1, P2, P3(s), P4, 'upper'); 
            end
            Ex_P1(s,k)  = (P3(s)./P2)* P4; % Expected number
            SetCounts(s,k) = P1;
            sNk = s+N*k;
            if mod(sNk,NK10th)==0 && sNk/NK10th < 11
                fprintf('%d%% ',  (sNk/NK10th) * 10);
            end
            
        end
        
    end
    fprintf('\n');
    
    
    %% FDR correction
    fprintf('Correction for Multiple Hypothesis testing...');
    Q_val = nan(size(P_val));
     if sum(~isnan(P_val(:))) < 2
        Q_val(:) = P_val(:);
    else
        NaNFiltp = (~isnan(P_val(:)));
        Q_val(NaNFiltp) = mafdr(P_val(NaNFiltp), 'BHFDR', true);
    end
    fprintf('done!\nFinalizing reports...');
    
    %% Reporting
    Input.GeneSets  = GeneSets;
    Input.GeneNames = GeneNames;
    Input.Clusters  = Clusters;
    Input.BackgroundSet = Clustered_Genes.Pej_ID;
    Results.P_val = P_val;
    Results.Q_val = Q_val;
    Results.GenesInClust        = GenesInClust;
    %
    try
        % Make an inverse Index for the clustered genes
        InvClustIdx = zeros(length(GeneNames),1);
        for iInv = 1:length(Clustered_Genes.Pej_ID)
            InvClustIdx(Clustered_Genes.Pej_ID(iInv)) = iInv;
        end
        avgEffectSize = nan([size(GenesInClust) size(Clustered_Genes.EffectSize,2)]);
        for s = 1:size(GenesInClust,1)
            for k = 1:size(GenesInClust,2)
                Xindx = InvClustIdx(GenesInClust{s,k});
                avgEffectSize(s,k,:)= mean(Clustered_Genes.EffectSize(Xindx,:));
            end
        end
        Results.meancoeffs = avgEffectSize;
        
    catch screwit
        disp('ups!')
    end
    Results.SetCountsInClusters = SetCounts;
    Results.SetCountsInBackground = P3;
    OverUnder = sign( Ex_P1 -  SetCounts);
    Results.OverUnderRepres = OverUnder;
    OV = {'Over' 'noth' 'Under'};
    
    Slash = find(Clustering_Output == '/', 1, 'last');
    if isempty(Slash); Slash=1;end;
    OutPutPrefix = [OutFldr '/' Clustering_Output(Slash+1:end) ' Test4EnrichmentDetails of '];
    OUT1 = [OutPutPrefix, ReferenceFilesLabels{SourceDB} '.mat'];
    save(OUT1, 'Results', 'Input' );
    fprintf('Done!\n')
    SigSets = sum((Q_val<=Thr),2)>0;
    
    MetaReport.DataBase(SourceDB,1) = ReferenceFilesLabels(SourceDB);
    MetaReport.TotalSetsInDB(SourceDB,1)            = size(Q_val,1);
    MetaReport.TotalTests(SourceDB,1)               = sum(~isnan(Q_val(:)));
    MetaReport.FDRcutOff(SourceDB,1)                = Thr;
    MetaReport.EnrichmentHits(SourceDB,1)           = sum(Q_val(:)<=Thr);
    MetaReport.Percent_H0_rejected(SourceDB,1)      = MetaReport.EnrichmentHits(SourceDB,1) ./sum(~isnan(Q_val(:)))*100;
    MetaReport.EnrichedSetsPerCluster(SourceDB,1)   = MetaReport.EnrichmentHits(SourceDB,1)/K;
    
    fprintf('There were %d significant sets under %d%% FDR\n', sum(SigSets),  Thr*100);
    disp(['with agerage of ' num2str(round(sum(sum((Q_val<=Thr)))/K)) ' sets per cluster']);
    fprintf('Results saved at: %s.\n\n', OUT1)
    %% make the overall report
    for s = 1:N
        if ~SigSets(s)
            continue
        end
        for k = 1:K
            if Q_val(s, k)<= Thr
                fprintf(F, '%s\t%s\t%s\t%s\t%2.1f\t%2.1f\t%d\t%.2f\t%d\t%.1f\t%s\t%d\t' , ...
                    ReferenceFilesLabels{SourceDB}, Clustered_Genes.ClusterLabels2{k}, GeneSets{s}.Name{1},OV{OverUnder(s,k)+2},...
                    log10(P_val(s,k)), log10(Q_val(s,k)), P3(s), Ex_P1(s,k),  SetCounts(s,k), Results.meancoeffs(s,k) , GeneSets{s}.Source{1}, k);
                %fprintf(F, '\t%s', GeneNames{intersect(currClust, GeneSets(s).Genes)});
                fprintf(F, '%s;', GeneNames{GenesInClust{s,k}});
                fprintf(F, '\n');
            end
        end
        
    end
end
fclose(F);
fprintf('Overall Report saved at: %s.\n', OverallReport)
Pej_Write_Table([OverallReport(1:end-4) '-MetaReport.txt'], MetaReport);
Ignore = system(['wc -l ' OverallReport]);
%PStatistics
end



function  Clustered_Genes = Get_PejIDs(Clustered_Genes, GeneNamesReference, Xref)
GeneIDs_new = Pej_Xref(Clustered_Genes.IDs, Xref);
[C, ai, bi] = intersect(GeneIDs_new.Associated_Gene_Name, GeneNamesReference);

tmp= Clustered_Genes.ClusterLabels2;
Clustered_Genes = rmfield(Clustered_Genes, 'ClusterLabels2');

Clustered_Genes = Pej_Struct_RowSelect(Clustered_Genes, ai);
Clustered_Genes.Pej_ID = bi;
[~, I]= sort(Clustered_Genes.Pej_ID);
Clustered_Genes = Pej_Struct_RowSelect(Clustered_Genes, I);

Clustered_Genes.ClusterLabels2 = tmp;
end