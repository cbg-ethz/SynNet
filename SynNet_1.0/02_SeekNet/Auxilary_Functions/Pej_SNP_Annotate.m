% SNP needs to have at least two fields, "position" and "chr".
% This implementation is very slow and stupid because I was lazy. you can
% just compare the start and end positions of the annotation to the
% position of SNP instead of using the more general genomic intersection
% script.

% Pejman April 2015, This code is incomplete! It finds the overlaps but
% then nothing
%-------------------

function SNP = Pej_SNP_Annotate(SNP, GTFAnnotationFile)

disp('Reading in annotation file...')
fid = fopen(GTFAnnotationFile, 'r');
InBuff = textscan(fid, '%s\t%s\t%*s\t%f\t%f\t%*s\t%*c\t%*s\t%[^\n]');fclose(fid);

Annotation.GeneID = Get_GTF_Atribute(InBuff{5}, 'gene_id');InBuff{5}=[];
Annotation.chr      = InBuff{1};InBuff{1}=[];
Annotation.GeneType = InBuff{2};InBuff{2}=[];
Annotation.start    = InBuff{3};InBuff{3}=[];
Annotation.end      = InBuff{4};InBuff{4}=[];
disp('done!')

Annotation.dir = repmat('+', length(Annotation.start),1);
SNP.dir        = repmat('+', length(SNP.position),1);
SNP.start      = SNP.position;
SNP.end        = SNP.position;

% If chr IDs are numeric make them cell strings
if isnumeric(SNP.chr)
    for i = length(SNP.chr):-1:1
        tmp{i,1} = num2str(SNP.chr(i));
    end
    SNP.chr = tmp;
end

[I]= Pej_Genomic_Intersection(SNP, Annotation);
SNP = rmfield(SNP, {'dir', 'start', 'end'});

T = Annotation.GeneType(I(:,2));
keepF = strcmp(T, 'protein_coding') | strcmp(T, 'lincRNA');
I = I(keepF,:);
PercUniq = 100*length(unique(I(:,1)))/length(SNP.position);
disp([num2str(PercUniq) '% were annotated'])

if length(unique(I(:,1)))<length((I(:,1)))
    warning('Some SNPS have more than one annotation, I kept the last annotation only!');
end
SNP.GeneID = cell(size(SNP.position));
SNP.GeneID(I(:,1)) = Annotation.GeneID(I(:,2));
end

function Value = Get_GTF_Atribute(S, Attr)
% Example:
% fid = fopen('Myfile.gtf', 'r');
% Annotatoin = Read_GTF(fid)
% Get_GTF_Atribute(Annotatoin.HitGeneType_LociType_Details(:,3), 'gene_name');
Value = regexp(S, [Attr ' "([^"]*)";'], 'tokens', 'once');
Value = vertcat(Value{:});
end