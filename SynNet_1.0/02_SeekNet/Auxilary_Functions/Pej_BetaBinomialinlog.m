function lpx = Pej_BetaBinomialinlog(x, xc, p, vScale)
vScale = exp(vScale);
if p==1
    lpx = log(double(xc==0));
    return
end

if p==0
    lpx = log(double(x==0));
    return
end

if vScale == Inf
    % binomial
    lpx = log(binopdf(x, x+xc, p));
    if ~isfinite(lpx)
        warning('Numerical failiure!')
    end
    return
end

if length(p)>1
    warning('this function is not written for vectors!')
end
a = vScale .* p;
b = a .* (1-p)./p;
lpx = (gammaln(x+xc + 1)-gammaln(x + 1)-gammaln(xc + 1)+...
    betaln((a + x),(b + xc))-betaln(a,b));

NaNflt = isnan(lpx);
if any(NaNflt)
    warning('Numerical failiure, substituted with binomial pdf!')
    lpx(NaNflt) = log(binopdf(x(NaNflt), x(NaNflt)+xc(NaNflt), p));
end
end