function [I1, I2]= Pej_Genomic_Intersection_2(Loci1, Loci2)
% This function reports the overlaps between two sets of genomic locations.

% NOTE: this function does THE SAME JOB as "Pej_Genomic_Intersection.m", it
% just has a different implementation. The other one is usually faster than
% this, but this one can be faster if there are a large amount of
% overlapping loci at one given point, like 100 or 1000 overlapping things at one point. 


% INPUTS:
% Two Loci structures:
% Loci.chr[nx1]     : chromosome names: (string)
% Loci.start[nx1]   : start of the position: (positive int)
% Loci.end[nx1]     : end   of the position: (positive int)
% Loci.dir[nx1]     : strand of the locale: (char: [ '+' / '-'])

% OUTPUTS:
% I1, I2, and two cell arrays such that: n-th locale in the Loci1 intersects
% with all locales in Loci2 with the index I1{n}


% Writtn by Pejman Mohammadi
% 26 June 2014, Basel
% Pejman.m@gmail.com
% ---------------------------------

L1 = size(Loci1.start,1);
L2 = size(Loci2.start,1);

if ~isfield(Loci1, 'chr_num')
    tmpCHRs  =[Loci1.chr; Loci2.chr];
    [~, ~, tmpCi] = unique(tmpCHRs);
    Loci1.chr_num = tmpCi(1:L1);
    Loci2.chr_num = tmpCi(L1+1:end);
end

Loci1.ID(:,1) = 1:L1;
Loci2.ID(:,1) = 1:L2;

senseFilt1 =  Loci1.dir == '+';
senseFilt2 =  Loci2.dir == '+';

I1 = cell(size(senseFilt1));
I2 = cell(size(senseFilt2));

fprintf('\n%d loci to compare against %d loci\n', L1, L2);
tic
fprintf('\nAnalyzing sense strand...\nChr: ');
[I1( senseFilt1), I2( senseFilt2) ] = GGI_noStrand(Pej_Struct_RowSelect(Loci1, senseFilt1), Pej_Struct_RowSelect(Loci2, senseFilt2));
fprintf('\nAnalyzing antisense strand...\nChr: ');
nS1 = ~senseFilt1;nS2 = ~senseFilt2;
[I1(nS1), I2(nS2) ] = GGI_noStrand(Pej_Struct_RowSelect(Loci1,nS1), Pej_Struct_RowSelect(Loci2,nS2));
fprintf('\nDone and Done! (%.0f secs)\n', toc);


end

function [I1, I2]= GGI_noStrand(Loci1, Loci2)
% This fuction works for inputs that come the same strand.
s1 = size(Loci1.chr_num);
s2 = size(Loci2.chr_num);

I1 = cell(s1);
I2 = cell(s2);


CHRs = unique(union(unique(Loci1.chr_num), unique(Loci2.chr_num))); % find all Chromosomes

unused1 = true(s1);
unused2 = true(s2);

for c = 1: length(CHRs)
    fprintf('%d ', CHRs(c));
    F1 = false(s1);
    F2 = false(s2);
    
    F1(unused1) = Loci1.chr_num(unused1) == CHRs(c);
    F2(unused2) = Loci2.chr_num(unused2) == CHRs(c);
    
    [I1(F1), I2(F2)] = GGI_noStrand_noChr(Pej_Struct_RowSelect(Loci1, F1), Pej_Struct_RowSelect(Loci2, F2));
    
    unused1(F1) = false;
    unused2(F2) = false;
end

end


function [I1, I2] = GGI_noStrand_noChr(Loci1, Loci2)
% This fuction works for inputs that come the same strand of the same chromosome.

if isempty(Loci1.start) || isempty(Loci2.start)
    I1=repmat({[]}, size(Loci1.start));
    I2=repmat({[]}, size(Loci2.start));
    return
end

if issorted(Loci1.start)
    if Loci1.start(1)<=Loci1.start(end)
        % it's ascending
        ix1 = 1:length(Loci1.start);
    else
        % it's descending
        ix1 = length(Loci1.start):-1:1;
    end
else
    [~, ix1] = sort(Loci1.start, 'ascend');
end

if issorted(Loci2.start)
    if Loci2.start(1)<=Loci2.start(end)
        % it's ascending
        ix2 = 1:length(Loci2.start);
    else
        % it's descending
        ix2 = length(Loci2.start):-1:1;
    end
else
    [~, ix2] = sort(Loci2.start, 'ascend');
end


[I1(ix1), I2(ix2)] = GGI_noStrand_noChr_sorted(Pej_Struct_RowSelect(Loci1, ix1), Pej_Struct_RowSelect(Loci2, ix2));
end

function [I1, I2] = GGI_noStrand_noChr_sorted(Loci1, Loci2)
% try
    % This fuction works for inputs that come the same strand of the same chromosome, and are sorted ascendingly for start positions.
L1 = size(Loci1.start,1);
L2 = size(Loci2.start,1);

I1 = cell(L1,1);
I2 = cell(L2,1);

mx = min([Loci1.start(1), Loci2.start(1)]);
Mx = max([max(Loci1.end), max(Loci2.end)]);

pointer1 = 1;
pointer2 = 1;

Overlapping1 = []; New_overlaps1 = [];
Overlapping2 = []; New_overlaps2 = [];

% these are the points that matter!
Milestones = union(union(Loci1.start, Loci2.start), union(Loci1.end, Loci2.end))';

for t = 1:length(Milestones)
    x = Milestones(t);
    % Check if overlapping hits are still overlapping
    Retain_OL1 = Loci1.end(Overlapping1)>=x;
    Retain_OL2 = Loci2.end(Overlapping2)>=x;
    
    % Check if new overlapping sets have started
    while pointer1<=L1 && x == Loci1.start(pointer1)
        New_overlaps1 = union(New_overlaps1, pointer1);
        pointer1 = pointer1+1;
    end
    
    while pointer2<=L2 && x == Loci2.start(pointer2)
        New_overlaps2 = union(New_overlaps2, pointer2);
        pointer2 = pointer2+1;
    end
    
    % Update Overlapping sets
    Overlapping1 = [Overlapping1(Retain_OL1), New_overlaps1]; New_overlaps1 = [];
    Overlapping2 = [Overlapping2(Retain_OL2), New_overlaps2]; New_overlaps2 = [];
    
    % Report Overlaps
    if ~isempty(Overlapping1) && ~isempty(Overlapping2)
        % Add the current overlapping loci to those that already exist
        for i = 1:length(Overlapping1)
            I1{Overlapping1(i)} = union(I1{Overlapping1(i)}, Overlapping2);
        end
        
        for i = 1:length(Overlapping2)
            I2{Overlapping2(i)} = union(I2{Overlapping2(i)}, Overlapping1);
        end
    end
    
    
end

% Translate the indexes to the original IDs
for i = 1:L1; I1{i} = Loci2.ID(I1{i}); end
for i = 1:L2; I2{i} = Loci1.ID(I2{i}); end    
% catch err
%     3
% end
end