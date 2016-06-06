function x = Pej_BetaBinomial_rnd(N, p, vScale)
if size(N,2)>1
    error('N is assumed to be a nx1 column')
end
vScale = exp(vScale);

a = vScale .* p;
b = a .* (1-p)./p;

n = length(N);
q = betarnd(a,b, n,1);
x = binornd(N,q, n,1 );

end
