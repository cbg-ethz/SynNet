% This functions makes a random circuit.
% MaxAnd, and Max or correspond to the size of the circiut, #Rows, and
% #Columns in the output matrix respectively.
% Ngene is the total number of available genes in the system.
function F = Build_Random_Function(Ngene, Consts)
addpath('Auxilary_Functions/');
if nargin<2; Consts = Fetch_Constraints(); end;
MaxAnd          = Consts.MaxAnd;
MaxOr           = Consts.MaxOr;
MaxOrCount      = Consts.MaxOrCount;
MaxNotCount     = Consts.MaxNotCount;
NotRate         = Consts.Simulation_NotRate;
MaxNG = min([
    MaxAnd*MaxOr
    Consts.MaxCircuitSize
    min(MaxOrCount+MaxNotCount, MaxAnd)+ min(MaxOrCount, MaxAnd) * (MaxOr-1)
    ]); % Maximum possible number of genes

if nargin==0; Ngene = MaxNG;end

NG = randi(MaxNG, 1); % Number of genes to be used (if a gene is used in two positions it's considered twice)
GI = [randi(Ngene,NG, 1); zeros(MaxAnd*MaxOr-NG,1)]; % the genes involved and the rest zero

Practical = false;C = 0;tmpNprn=0;
fprintf('Synthetizing a valid circuit with %d genes... (Try no: ', NG);
while ~(Practical)
    fprintf(repmat('\b', 1, tmpNprn));
    tmpNprn = fprintf('%d)',C);
    C = C+1;
    F = Shape_Circuit(GI, MaxAnd, MaxOr, NotRate);
    Practical = Ispractical(F,'Discard');
end
F = Refine_HitPool(F);
Print_Function(F);

end

function F = Shape_Circuit(GI, MaxAnd, MaxOr, NotRate)

F = reshape(GI(randperm(MaxAnd*MaxOr)),MaxAnd, MaxOr);  % put them st random in a function
% the next line puts negative sign in front of half of the single gene ORs(rows)
F = F .* (ones(size(F)) - 2 * ones(size(F)) .* repmat(rand(MaxAnd, 1)>NotRate, 1, MaxOr) .* repmat(sum(F>0,2)==1, 1, MaxOr));
end