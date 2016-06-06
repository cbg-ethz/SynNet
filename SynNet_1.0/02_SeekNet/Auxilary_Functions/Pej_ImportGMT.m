% This is the version 3. Functionality is exactly the same as before, it's
% just faster  :)
% Pejman 26 Feb 2013 Lausanne, CHUV
%---------
% Reading in GSEA MSigDB files in .gmt format
% this file gets a .gmt file available here:
% http://www.broadinstitute.org/gsea/msigdb/collections.jsp
% and a list of Gene names, then it gives out an array of genesets describing who is
% connected to who, using the index of gene in the GeneNames input
% parameter.

function Pej_ImportGMT(AllGenes, DB_Path)

if nargin < 2
    DB_Path = '/Users/pejmanm/Desktop/LocalTMP/PEJ_Resources/GeneSets/';
end

if nargin < 1
    tmpAllGenes = Pej_Read_Table('/Users/pejmanm/Desktop/LocalTMP/PEJ_Resources/Gene-Xref-25Apr2015.txt', [], false);
    AllGenes.GeneIDs = unique(tmpAllGenes.Associated_Gene_Name, 'sorted'); clear tmpAllGenes
    AllGenes.Pej_ID  = 1:length(AllGenes.GeneIDs);
end
GeneNames = AllGenes.GeneIDs;
% build an index for the first Char of the gene
[GN GI]  = sort((GeneNames));
for gi = length(GI):-1:1
    tmpGNindex(gi) = GN{gi}(1);
end

N = length(GeneNames);

Pathbkp = pwd;
cd(DB_Path);
listing = dir();
for dbFile = 1:length(listing)
    [Fpath FName Fext] = fileparts(listing(dbFile).name);
    if ~isempty(FName) && FName(1)~= '.' && strcmpi(Fext, '.gmt')
        %% It's a .gmt file!
        T = tic;
        disp(['Importing ' listing(dbFile).name '...'])
        FilePath = listing(dbFile).name;
        if exist([FilePath ' GeneSets.mat'], 'file')
            %% It's there already
            Prev = load([FilePath ' GeneSets.mat'], 'GeneNames');
            if isequaln(Prev.GeneNames, GeneNames)
                disp('Import Skipped: Index Already existed!')
                continue;
            end
        end
        
        %% Import the file
        %get the length of the file
        FilePathF = which(FilePath);
        FilePathF = strrep(FilePathF, ' ', '\ ');
        [~, numlines] = system( ['wc -l ', FilePathF] );
        numlines = sscanf(numlines, '%f ');
        disp([num2str(numlines) ' sets found.'])
        numlinesperc = round(numlines/100);
        F = fopen(FilePath, 'r');
        Lcnt = 0;tic
        percrep=0;
        for Lcnt = numlines:-1:1
            L = fgetl(F);
            tmpSet = regexp(L,'\t','split');
            GeneSet.Name   = tmpSet(1);
            GeneSet.Source = tmpSet(2);
            GeneSet.SourceFile = FilePath;
            tmpSet = unique(tmpSet(3:end), 'sorted');
            %             [~,~,tmpSetib] = intersect(tmpSet,GN); % (A,B)
            tmpSetib = Pej_Intersect_SortedVectors(tmpSet,GN); % (A,B)
            GeneSet.Genes  = GI(tmpSetib)';
            GeneSets{numlines - Lcnt + 1} = GeneSet;
            clear GeneSet
            if floor(Lcnt / numlinesperc) == Lcnt / numlinesperc
                fprintf(repmat('\b', 1,percrep));
                percrep = fprintf('%d%% ',Lcnt/numlinesperc);
            end
        end
        fclose(F);
        for GS = 1:length(GeneSets)
            GeneSets{GS}.Genes_PejIDs = AllGenes.Pej_ID(GeneSets{GS}.Genes);
        end
        save([FilePath ' GeneSets.mat'], 'GeneSets', 'GeneNames');
        clear GeneSets
        
        T = toc(T);
        disp(['done! (' num2str(T) 'secs)'])
        
    end
end
cd(Pathbkp);
end