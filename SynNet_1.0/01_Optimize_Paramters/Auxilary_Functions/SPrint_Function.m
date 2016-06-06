function FunctionString = SPrint_Function(Function, Labels, Performance, FileID)
if nargin<4
    FileID = 1;
end
if nargin<2 || isempty(Labels)
    %% Make up Gene identifiers.
    for i = max(abs(Function(:))):-1:1
        Labels{i} = ['G_{' int2str(i) '}'];
    end
end

OutBuff{1}='';
Function = squeeze(Function);
if length(size(Function)) == 3
    FunctionString ='';
    for k = 1:min(Inf, size(Function,1))
        OutBuff{end+1}=sprintf('\nHit No. %d\n', k);
        if nargin>2
            OutBuff{end+1} = SPrint_Function(squeeze(Function(k,:,:)), Labels, Performance(k), FileID);
        else
            OutBuff{end+1} = SPrint_Function(squeeze(Function(k,:,:)), Labels);
        end
    end
else
    
    
    MaxAnd = size(Function,1);
    MaxOr  = size(Function,2);
    
    
  %  OutBuff{end+1}=sprintf('\n=============================\n');
    if nargin>2
        OutBuff{end+1}=sprintf('Performance: %.2f%%\n', Performance);
    end
    FirstLine = true;
    for i = 1:MaxAnd
        if ~all(Function(i,:)==0)
            if ~FirstLine;   OutBuff{end+1}=sprintf('AND ');     end;
            Newline = true;
            for j = 1:MaxOr
                tmpGI = Function(i,j);
                if tmpGI ~=0
                    if ~Newline;      OutBuff{end+1}=sprintf('|');      end
                    if tmpGI < 0;     OutBuff{end+1}=sprintf('n');        end
                    OutBuff{end+1}=sprintf('%s', Labels{abs(tmpGI)});
                    Newline = false;
                end
            end
            OutBuff{end+1}=sprintf(' ');
            FirstLine = false;
        end
    end
    OutBuff{end+1}=sprintf('\n');
end
FunctionString =sprintf('%s', OutBuff{:});
end