% This file gets a structure containing names and expressison for miRNAs,
% and a path to file containing mature miRNA sequences. It uses the
% sequence similarity to add up expression of similar miRNAs to each other.
% Additionally this code picks only the first gene from a list of
% equivalent ones and sets all the rest expression to zero in all samples.

% NOTE: The expressions should be Normalized before this file, an should
% not be renormalized afterwards.
%Pejman 2013 Jan
%------------------------------

function    miRNA_Data = Combine_Similar_miRNAs(miRNA_Data, Consts)
E = miRNA_Data.Dnormed;
L = lower(miRNA_Data.GeneNames);
N = length(L);
disp('Loading mature RNA sequences to Check similarities')
Seqs = fastaread(Consts.Mature_miRNAs, 'TrimHeaders' , true);
miLs = lower({Seqs(:).Header}');
miSq = {Seqs(:).Sequence}'; clear Seqs

[C ia ib] = intersect(L, miLs);

Ltmp = regexp(L, '/', 'split', 'once');
for l = 1:N
    if length(Ltmp{l})>1 && isempty(find(l==ia))
        % the   ID is a combined ID like: 'hsa-let-7a/hsa-let-7b'
        [Cl, ~, ibl] = intersect(Ltmp{l}, miLs);
        if ~isempty(Cl)
            C(end+1)  = L(l);
            ia(end+1) = l;
            ib(end+1) = ibl(1);
        end
    end
end

M =length(C);
disp([int2str(M) ' out of ' int2str(N) ' miRNA IDs were matched to the provided mature miRNAs'])
miRNA_Data.Mature_Sequence = cell(N,1);
miRNA_Data.Mature_Sequence(ia) = miSq(ib);

for i = M:-1:1
    Seeds(i,:) = miRNA_Data.Mature_Sequence{ia(i)}(2:8);
end

dSeeds = double(Seeds);

tmpEia = E(ia,:);
Eia    = nan(size(tmpEia));
Exclude = false(size(Eia,1),1);
% Add expressions of similar miRNAs together
for i = M:-1:1
    tmpD = sum(dSeeds ~= repmat(dSeeds(i,:),M,1),2); % Hamming distance to all seeds
    tmpF = tmpD<=Consts.Max_Seed_Similatiry;
    Eia(i,:) = sum(tmpEia(tmpF, :),1);
    SourceMiRia(:,i) =  tmpF;
    tmpS = sprintf('%s\t=', L{ia(i)});
    OutBuff = [L(ia(tmpF)), cellstr(Seeds(tmpF,:))]';
    GNia{i,1} = sprintf('%s+', L{ia(tmpF)});
    GNia{i}=GNia{i}(1:end-1); % remove the extra "+" from the end
    
    tmpS = [tmpS sprintf('\t%s(%s)', OutBuff{:,:})];
    %     if sumNaN(tmpEia(i,:))<(.5*sumNaN(Eia(i,:)))
    if any(sumNaN(tmpEia(i,:),2)<sumNaN(tmpEia(tmpF,:),2))
        % Exclude the miRNA from the analysis. Too low expression share among the neighbouring similar miRNAs
        Exclude(i)= true;
        tmpS = [tmpS sprintf('\t(Minor share - Excluded)')];
    else
        % Keep it!
    end
    FoutB{i} = [tmpS sprintf('\n')];
end

Eia(Exclude,:)=0; % Put all of them to zero expression
fprintf('%d miRNAs were excluded, due to low expression share among the neighbouring similar miRNAs.\n', sum(Exclude));

Fout = fopen([Consts.OutFolder '/Combine_Similar_miRNAs.log'], 'w');
fprintf(Fout, '%s', FoutB{:});
E(ia,:) = Eia; clear Eia
miRNA_Data.GeneGroupNames = cell(N,1);
miRNA_Data.GeneGroupNames(ia)=GNia; clear GNia
SourceMiR(:,ia) = SourceMiRia;  clear SourceMiRia

%% Remove redundancies
RdnFilt = false(M,1);
for i = 1:M-1
    for j = (i+1):M
        if isequaln(SourceMiR(:,i), SourceMiR(:,j))
            RdnFilt(j)= true;
            break
        end
    end
end

fprintf('%d miRNAs were redundant.\n', sum(RdnFilt));
% Remove the redundant ones
E(RdnFilt,:)=0; % Put all the redundant gene to zero expression

miRNA_Data.Dnormed = E;
fclose(Fout);
end

function S = sumNaN(X,d)
if nargin==1
    if size(X,1)>1
        d=1;
    else
        d=2;
    end
end
X(isnan(X))=0;
S = sum(X,d);
end
