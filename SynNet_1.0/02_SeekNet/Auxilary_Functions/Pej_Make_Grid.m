% This function gets a set of Domain Limits of size(n*2) like:
% [ 0 3; 2 5; 1 4];
% defining lower and uppper bound for each of the n dimentions.
% and produces a regular grid of "Density" in "Grid", with dimentions
% (Density.^n, n);
% May 2014 Pejman, pejman.m@gmail.com
%----------------------

function Grid = Pej_Make_Grid(DomainLimits, Density)

PDim = size(DomainLimits,1);
Grid = nan(Density .^ PDim, PDim);
for p = 1: PDim
    tmpP = repmat(linspace(DomainLimits(p,1), DomainLimits(p,2), Density), Density .^ (PDim-p),1);
    Grid(:,p) = repmat(tmpP(:), Density .^ (p-1),1);
end

end