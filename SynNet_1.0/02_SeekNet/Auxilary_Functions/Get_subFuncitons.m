% Warning: This Code is written assuming that the function would be small
% (~10 genes), it produces exponential amount of output!

function Sub_F = Get_subFuncitons(Function)
Function = squeeze(Function);
GL = find(Function~=0);
N = length(GL);
if N == 0
    Sub_F = [];
end

if N>10
    warning('Pruning function is written to perform an exhaustive search, this costs of O(2^Number of Genes in the circuit)!')
end
Exclude = Pej_Make_Grid(repmat([0,1],N,1), 2)== 0;


for i = size(Exclude,1):-1:1
    tmpF = Function;
    tmpF(GL(Exclude(i,:)))=0;
    Sub_F(i,:,:) = tmpF;
end
end