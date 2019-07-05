function xst = softth(x,l,nonneg)
% soft threshold of x by l;
if nargin==2;
    nonneg=0;
end
xst=sign(x).*max(0, abs(x)-l);
if nonneg==1;
    xst(xst<0) = 0;
end