function [a, d, b] = ordinalcomponentupdate(x,a,b,D,l,nonneg)

maxit = 5; 
conv =1e-6;
Fit = 0;
diff = 1;
c = 0;
SSETOT = norm(x,'fro');
while diff > conv & c<maxit
    c = c+1;
    % perform ordinal clustering over samples based on least squares ordering
    [aa, id] = sort(a);
    [~,id2] = sort(id);
    da = aa(2:end) - aa(1:end-1);
    da = 1./(da + 0.1);
    da = da / mean(da);
    % update a
    a = genSMRupdate_dual(x(id,:)',b,diag(da)*D,l);
    a = colNorm(a(id2));
    
    % update b by LS
    b = colNorm(x'*a);
    b = colNorm(searchL_BIC(x,a,nonneg)); 
    d = a'*x*b;
    
    % check fit
    Fitold = Fit;
    SSE = norm(x - d*a*b','fro');
    Fit = (SSETOT - SSE)/SSETOT;
    diff = (Fit - Fitold);
    FIT(c) = Fit;
    sse(c) = SSE;
end

