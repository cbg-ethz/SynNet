function Print_Function(Function, Labels, Performance, FileID)
if nargin<4
    FileID = 1;
end
if nargin<2 || isempty(Labels)
    %% Make up Gene identifiers.
    for i = max(abs(Function(:))):-1:1
        Labels{i} = ['G' int2str(i)];
    end
end

Function = squeeze(Function);
if length(size(Function)) == 3
    for k = 1:min(5, size(Function,1))
        fprintf(FileID,'\nHit No. %d', k);
        if nargin>2
            Print_Function(squeeze(Function(k,:,:)), Labels, Performance(k), FileID)
        else
            Print_Function(squeeze(Function(k,:,:)), Labels)
        end
    end
else
    
    
    MaxAnd = size(Function,1);
    MaxOr  = size(Function,2);
    
    
    fprintf(FileID,'\n=============================\n');
    if nargin>2
        fprintf(FileID,'Performance: %.2f%%\n', Performance);
    end
    FirstLine = true;
    for i = 1:MaxAnd
        if ~all(Function(i,:)==0)
            if ~FirstLine;   fprintf(FileID,'-AND-\n');     end;
            Newline = true;
            for j = 1:MaxOr
                tmpGI = Function(i,j);
                if tmpGI ~=0
                    if ~Newline;      fprintf(FileID,' OR ');      end
                    if tmpGI < 0;     fprintf(FileID,'NOT');        end
                    fprintf(FileID,'(%s)', Labels{abs(tmpGI)});
                    Newline = false;
                end
            end
            fprintf(FileID,'\n');
            FirstLine = false;
        end
    end
    fprintf(FileID,'\n');
end
end