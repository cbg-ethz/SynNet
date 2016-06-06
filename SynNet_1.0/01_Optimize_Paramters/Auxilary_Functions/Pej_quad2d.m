function X = Pej_quad2d(dx, D1L, D1H, D2L, D2H, RelativeTol, BreakPoint)
SymmetricDomain = all(BreakPoint == mean([D1L, D1H; D2L, D2H],2)); % Matlab's multidimenstional integrators work very fast on symmetric domains

if SymmetricDomain
    X = quad2d(dx, D1L, D1H, D2L, D2H,'AbsTol', 0, 'RelTol', RelativeTol);
else
    BrP = BreakPoint;
    EPS = 0.001;
    X = ...
        quad2d(dx, D1L      , BrP, D2L      , BrP,'AbsTol', 0, 'RelTol', RelativeTol) + ...
        quad2d(dx, BrP + EPS, D1H, D2L      , BrP,'AbsTol', 0, 'RelTol', RelativeTol) + ...
        quad2d(dx, D1L      , BrP, BrP + EPS, D2H,'AbsTol', 0, 'RelTol', RelativeTol) + ...
        quad2d(dx, BrP + EPS, D1H, BrP + EPS, D2H,'AbsTol', 0, 'RelTol', RelativeTol);
end