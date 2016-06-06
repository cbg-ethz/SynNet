function P = BackProb_log(X)
%% UN-NORMALIZED Background probability distribution function of the LOG-transformed input
for d = size(X,1):-1:1; Pd(d,:) = MarginalProb_log(X(d,:), d);end
P = prod(Pd,1);

end
function P = MarginalProb_log(x,d)
%% Marginal prior probability distribution  function of log10(x) for dimension "d"
global PDname Pa Pb Pc Pscalefactor

P = pdf(PDname{d},x,Pa(d), Pb(d), Pc(d)) * Pscalefactor(d) ;
end