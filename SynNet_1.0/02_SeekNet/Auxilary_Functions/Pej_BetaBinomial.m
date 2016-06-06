function px = Pej_BetaBinomial(x, xc, p, vScale)
vScale = exp(vScale);

a = vScale * p;
b = a * (1-p)/p;
px = exp(gammaln(x+xc + 1)-gammaln(x + 1)-gammaln(xc + 1)+...
    betaln((a + x),(b + xc))-betaln(a,b));
end