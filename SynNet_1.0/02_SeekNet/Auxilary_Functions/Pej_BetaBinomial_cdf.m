function px = Pej_BetaBinomial_cdf(x, xc, p, vScale)
%vScale = exp(vScale);
%a = vScale * p;
%b = a * (1-p)/p;

px = zeros(size(x));
for k = 1:length(x)
    px(k) = Pej_BetaBinomial_cdf_single(x(k), xc(k), p, vScale);
end
end


function px = Pej_BetaBinomial_cdf_single(x, xc, p, vScale)
% We assume X < Xc always
if x>xc
    t=xc;
    xc=x;
    x=t;
end

px = 0;
% px2 = 0;

vScale2 = exp(vScale);
a = vScale2 * p;
b = a * (1-p)/p;
t1 = gammaln(x+xc + 1);
t2 = betaln(a,b);
for i = 0:x
    txc = xc+x-i;
%     px = Pej_BetaBinomial(i, xc+x-i, p, vScale) +px; 
    px = exp(t1-gammaln(i + 1)-gammaln(txc + 1)+...
    betaln((a + i),(b + txc))-t2) + px;
end
if px>.5
    px = 1-px;
end
px = px*2;% Two tailed


if px>1
    1
end
end