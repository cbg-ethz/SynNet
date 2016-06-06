function X = Pej_Normalize_Columns(X,Offset )
XL = log(X);     clear X
if nargin<2
    ref = mean(XL,2);
    Xd = XL - repmat(ref,1, size(XL,2));
    Filt = ~any(isnan(Xd),2);
    Offset = median(Xd(Filt,:));
end
X = exp(XL - repmat(Offset, size(XL,1),1));
end