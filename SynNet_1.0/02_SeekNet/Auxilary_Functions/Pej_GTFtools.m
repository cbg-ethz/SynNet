%% split files
% First I broke the files into Chromosom-strands with the followoing two
% lines; Ideally this two lines need to have one "if" in them to close the
% files when not needed otherwise you might get way too many files open at
% the same time
% awk '{print >"CH_"$1$7".txt"}' Homo_sapiens.GRCh37.72.gtf


function Annotatoin = Read_GTF(fid)
InBuff = textscan(fid, '%*s\t%s\t%s\t%d\t%d\t%*s\t%*c\t%*s\t%[^\n]');fclose(fid);

GeneType = InBuff{1};
LociType = InBuff{2};
Details  = InBuff{5};
Annotatoin.LociCord = [InBuff{3} InBuff{4}]; clear InBuff
Annotatoin.GeneType_LociType_Details = [GeneType LociType Details];

end

function Value = Get_GTF_Atribute(S, Attr)
% Example:
% fid = fopen('Myfile.gtf', 'r');
% Annotatoin = Read_GTF(fid)
% Get_GTF_Atribute(Annotatoin.HitGeneType_LociType_Details(:,3), 'gene_name');
Value = regexp(S, [Attr ' "([^"]*)";'], 'tokens', 'once');
Value = vertcat(Value{:});
end



