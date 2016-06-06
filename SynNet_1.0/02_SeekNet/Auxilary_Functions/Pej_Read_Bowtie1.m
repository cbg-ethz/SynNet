% This just reads Default bowtie output:
% bowtie outputs one alignment per line. Each line is a collection of 8 fields separated by tabs; from left to right, the fields are:
%
% 1 Name of read that aligned.
% Note that the [SAM specification] disallows whitespace in the read name. If the read name contains any whitespace characters, Bowtie 2 will truncate the name at the first whitespace character. This is similar to the behavior of other tools.
%
% 2 Reference strand aligned to, + for forward strand, - for reverse
%
% 3 Name of reference sequence where alignment occurs, or numeric ID if no name was provided
%
% 4 0-based offset into the forward reference strand where leftmost character of the alignment occurs 
% NOTE >>>>I make the position 1-based! <<<<<

% 5 Read sequence (reverse-complemented if orientation is -).
% If the read was in colorspace, then the sequence shown in this column is the sequence of decoded nucleotides, not the original colors. See the Colorspace alignment section for details about decoding. To display colors instead, use the --col-cseq option.
%
% 6 ASCII-encoded read qualities (reversed if orientation is -). The encoded quality values are on the Phred scale and the encoding is ASCII-offset by 33 (ASCII char !).
% If the read was in colorspace, then the qualities shown in this column are the decoded qualities, not the original qualities. See the Colorspace alignment section for details about decoding. To display colors instead, use the --col-cqual option.
%
% 7 If -M was specified and the prescribed ceiling was exceeded for this read, this column contains the value of the ceiling, indicating that at least that many valid alignments were found in addition to the one reported.
% Otherwise, this column contains the number of other instances where the same sequence aligned against the same reference characters as were aligned against in the reported alignment. This is not the number of other places the read aligns with the same number of mismatches. The number in this column is generally not a good proxy for that number (e.g., the number in this column may be '0' while the number of other alignments with the same number of mismatches might be large).
%
% 8 Comma-separated list of mismatch descriptors. If there are no mismatches in the alignment, this field is empty. A single descriptor has the format offset:reference-base>read-base. The offset is expressed as a 0-based offset from the high-quality (5') end of the read.

% Reference: http://bowtie-bio.sourceforge.net/manual.shtml#sam


% Pejman 2014
%--------------

function Data = Pej_Read_Bowtie1(Path2File)

DLM = '[\t]'; % list of potential delimiters
fprintf(['Input file: %s'], Path2File);

FormatS = '%s%c%s%f%s%s%f%s';
Header = {'Seq_ID', 'dir', 'chr', 'start', 'Seq', 'Quals', 'Other_aligns', 'Mismatches'};
%% Read file
Fin = fopen(Path2File, 'r');
DES_Res = textscan(Fin, FormatS,'Delimiter', DLM, 'commentStyle','@');%, 'bufsize', 1E+5);
for i = 1:length(DES_Res)
    Data.(Header{i}) = DES_Res{i}; DES_Res{i} = [];
end
fclose(Fin);
Data.start = Data.start+1; % I make the position 1-based!
fprintf('\tdone!\n')
end



