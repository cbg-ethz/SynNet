% This functiobs joins structures based on a one or more Keys
function [InputArray] = Pej_Struct_Join(KeysFields, InputArray)

%% Make the Key format
KeyFormat = '';
for i = length(KeysFields):-1:1
    if isnumeric(InputArray{1}.(KeysFields{i})) || islogical(InputArray{1}.(KeysFields{i}))
        if isequalwithequalnans(InputArray{1}.(KeysFields{i}), round(InputArray{1}.(KeysFields{i})))
            %print as integer
            KeyFormat= ['%d_' KeyFormat ];
        else
            % print as float
            KeyFormat= ['%f_' KeyFormat ];
        end
        KeyType(i) = 1; % number
    else
        KeyFormat= [KeyFormat '%s_'];
        KeyType(i) = 2; % string
    end
end
KeyFormat = [KeyFormat '\n'];

%% build the key for each structure
for K = 1:length(InputArray)
    clear KeyBuff
    for i = length(KeysFields):-1:1
        if KeyType(i) == 1
            KeyBuff(:,i)  = num2cell(InputArray{K}.(KeysFields{i}));
        else
            KeyBuff(:,i)  = (InputArray{K}.(KeysFields{i}));
        end
    end
    KeyBuff   = KeyBuff';
    InputArray{K}.Pej_Key = regexp(sprintf(KeyFormat, KeyBuff{:}), '\n', 'split')';
    if isempty(InputArray{K}.Pej_Key{end})
        % regexp has added one extra cell at the end
        InputArray{K}.Pej_Key(end)=[];
    end
    
    if length(unique(InputArray{K}.Pej_Key))<length(InputArray{K}.Pej_Key)
        error('Your set of columns do not make a unique key!')
    end
    
    %% Make a list of all keys
    % Here I discard all the keys that are not observed in all samples
    if K ==1
        AllKeys = InputArray{K}.Pej_Key;
    else
        AllKeys = intersect(AllKeys, InputArray{K}.Pej_Key);
    end
end
%% Compare everything to the All Keys and return the lists in a comparable index
for K = 1:length(InputArray)
    [~, ~, ib] = intersect(AllKeys, InputArray{K}.Pej_Key);
    InputArray{K} = Pej_Struct_RowSelect(InputArray{K}, ib);
end


end
