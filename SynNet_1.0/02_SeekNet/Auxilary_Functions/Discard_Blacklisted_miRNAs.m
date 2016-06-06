function miRNA_Data = Discard_Blacklisted_miRNAs(miRNA_Data, Consts)

Fin = fopen(Consts.MiRNAs_Blacklist,'r');
if Fin == -1
    error(['The blacklist file address: ' Consts.MiRNAs_Blacklist ' is invalid.']);
end
Blklist = textscan(Fin, '%s');
Blklist = lower(Blklist{1});
L = lower(miRNA_Data.GeneNames);
[C,ia, ~] = intersect(L, Blklist);
miRNA_Data.Dnormed(ia,:)=0;
fprintf('Out of %d miRNAs in blacklist, %d were matched to the data and were excluded.\n', length(Blklist), length(C));
end