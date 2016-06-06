function s = Pej_Decode_NTs(s)
s = strrep(s,'1','A');
s = strrep(s,'2','C');
s = strrep(s,'3','G');
s = strrep(s,'4','T');
s = strrep(s,'5','N');
end