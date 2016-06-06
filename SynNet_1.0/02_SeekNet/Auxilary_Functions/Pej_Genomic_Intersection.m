function [I]= Pej_Genomic_Intersection(Loci1, Loci2)
% This function reports the overlaps between two sets of genomic locations.

% NOTE: this function does THE SAME JOB as "Pej_Genomic_Intersection_2.m", it
% just has a different implementation. The other one is usually slower than
% this, but this one can be slow if there are a large amount of
% overlapping loci at one given point, like 100 or 1000 overlapping things at one point. 


% INPUTS:
% Two Loci structures:
% Loci.chr[nx1]     : chromosome names: (string)
% Loci.start[nx1]   : start of the position: (positive int)
% Loci.end[nx1]     : end   of the position: (positive int)
% Loci.dir[nx1]     : strand of the locale: (char: [ '+' / '-'])
% Loci.chr_num[nx1] : [optional] chromosome numeric Ids: (number)

% OUTPUTS:
% I: A 2xz array such that each row, [n, m], says n-th locale in the Loci1 intersects
% with m-th locale in the loci 2. If one set has more than one intersecting
% sets in the other set, it will show up in multiple lines.


% Writtn by Pejman Mohammadi
% 27 June 2014, Basel
% Pejman.m@gmail.com
% ---------------------------------

% An ID field is defined in this code so we remove it from the data if it happens to be there!
if isfield(Loci1, 'ID');  
    warning('There was an ID field in the data, FYI I removed it :)')
    Loci1 = rmfield(Loci1, 'ID');
end
if isfield(Loci2, 'ID');    
    warning('There was an ID field in the data, FYI I removed it :)')
    Loci2 = rmfield(Loci2, 'ID');
end

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

fprintf('\n%d loci to compare against %d loci\n', L1, L2);
tic
fprintf('\nAnalyzing sense strand...\nChr: ');
I =[];
[tmpI] = GGI_noStrand(Pej_Struct_RowSelect(Loci1, senseFilt1), Pej_Struct_RowSelect(Loci2, senseFilt2));
I = [I; tmpI];

fprintf('\nAnalyzing antisense strand...\nChr: ');
nS1 = ~senseFilt1;nS2 = ~senseFilt2;
[tmpI ] = GGI_noStrand(Pej_Struct_RowSelect(Loci1,nS1), Pej_Struct_RowSelect(Loci2,nS2));
fprintf('\nDone and Done! (%.0f secs)\n', toc);
I = [I; tmpI];


end

function [I]= GGI_noStrand(Loci1, Loci2)
% This fuction works for inputs that come the same strand.
s1 = size(Loci1.chr_num);
s2 = size(Loci2.chr_num);

I = [];

CHRs = unique(union(unique(Loci1.chr_num), unique(Loci2.chr_num))); % find all Chromosomes

unused1 = true(s1);
unused2 = true(s2);

for c = 1: length(CHRs)
    fprintf('%d ', CHRs(c));
    F1 = false(s1);
    F2 = false(s2);
    
    F1(unused1) = Loci1.chr_num(unused1) == CHRs(c);
    F2(unused2) = Loci2.chr_num(unused2) == CHRs(c);
    
    [tmpI] = GGI_noStrand_noChr(Pej_Struct_RowSelect(Loci1, F1), Pej_Struct_RowSelect(Loci2, F2));
    I = [I; tmpI];
    
    unused1(F1) = false;
    unused2(F2) = false;
end

end


function [I] = GGI_noStrand_noChr(Loci1, Loci2)
% This fuction works for inputs that come the same strand of the same chromosome.
I = [];
if isempty(Loci1.start) || isempty(Loci2.start)
    
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


[tmpI] = GGI_noStrand_noChr_sorted(Pej_Struct_RowSelect(Loci1, ix1), Pej_Struct_RowSelect(Loci2, ix2));
I = [I; tmpI];
end

function [I] = GGI_noStrand_noChr_sorted(Loci1, Loci2)
try
    % This fuction works for inputs that come the same strand of the same chromosome, and are sorted ascendingly for start positions.
    L1 = size(Loci1.start,1);
    L2 = size(Loci2.start,1);
    mL = ceil(mean(L1,L2)/10);
    
    I = zeros(mL, 2);
    
    pointer1 = 1;
    pointer2 = 1;
    
    Overlapping1 = []; New_overlaps1 = [];
    Overlapping2 = []; New_overlaps2 = [];
    
    % these are the points that matter!
    Milestones = union(union(Loci1.start, Loci2.start), union(Loci1.end, Loci2.end))';
    lenI = 0;
    for t = 1:length(Milestones)
        x = Milestones(t);
        % Check if overlapping hits are still overlapping
        Overlapping1(Loci1.end(Overlapping1)<x)=[];
        Overlapping2(Loci2.end(Overlapping2)<x)=[];
        
        % Check if new overlapping sets have started
        while pointer1<=L1 && x == Loci1.start(pointer1)
            New_overlaps1 = [New_overlaps1; pointer1];
            pointer1 = pointer1+1;
        end
        
        while pointer2<=L2 && x == Loci2.start(pointer2)
            New_overlaps2 = [New_overlaps2; pointer2];
            pointer2 = pointer2+1;
        end
        
        % Update Overlapping sets
        if ~isempty(New_overlaps1)
            Overlapping1 = [Overlapping1; New_overlaps1]; 
            New_overlaps1 = [];
        end
        if ~isempty(New_overlaps2)
            Overlapping2 = [Overlapping2; New_overlaps2]; 
            New_overlaps2 = [];
        end
        
        % Report Overlaps
        if ~isempty(Overlapping1) && ~isempty(Overlapping2)
            tmpI = Cart_prod(Overlapping1, Overlapping2);
            tL   = size(tmpI,1);
            if (lenI + tL)>size(I,1);
                %increase the size of I; This is just to pre-allocate I to speed up
                I(end+mL,:)=0;
            end
            I(lenI+(1:tL),:) = tmpI;
            lenI = lenI+tL;
        end
        
    end
I = unique(I(1:lenI,:), 'rows') ;
I(:,1) =  Loci1.ID(I(:,1));
I(:,2) =  Loci2.ID(I(:,2));
catch
2    
end

end

function P = Cart_prod(S1, S2)
if isscalar(S1) && isscalar(S2); P = [S1, S2]; return; end
S1 = unique(S1);S2 = unique(S2);
[X, Y] = meshgrid(S1, S2);
P = [X(:), Y(:)];
end
