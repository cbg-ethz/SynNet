function [X, NormalizedSamples]= Pej_Normalize_Columns_Robust(X,Offset )

F = sum(X>1,2)>(size(X,2)/2);
XL = log(X(F,:));
NormalizedSamples = sum(isfinite(XL))>(.75*size(XL,1)); % Only include samples that have over 75% of the expressed genes expressed in them!
XL = XL(:,NormalizedSamples);
ref = median(XL,2);
Xd = XL - repmat(ref,1, size(XL,2));
Filt = ~any(isnan(Xd),2);
Offset = 1./exp(median(Xd(Filt,:)));

X(:, NormalizedSamples) = X(:,NormalizedSamples) .* repmat(Offset, size(X,1),1);
X(:,~NormalizedSamples) = nan;
end