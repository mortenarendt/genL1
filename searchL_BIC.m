function bs = searchL_BIC(X,a,nonneg)
% l is the lambda softthresholding penalty for shrikage of Ails. 
% X is data
% a is the left component
n = numel(X);
m = 30; 
% calculate the LS estimate of b 
b = colNorm(X'*a); 
L = linspace(0,max(abs(b)),m); 
for i=1:m; 
    % shrink the vector by L(i)
    bs = softth(b,L(i),nonneg); 
    d = a'*X*bs;
    % calculate current loss
    E = X - a*d*bs';
    sse = sum(E(:).^2);
    BIC2(i,1) = log(sse / n);  
    BIC2(i,2) = log(n)/n*sum(bs~=0); 
end
BIC = sum(BIC2'); 
% plot(L,BIC,'bo-',L,BIC2','--'); shg
% plot(L,BIC,'bo-'); shg
[mBIC id] = min(BIC); 
bs = softth(b,L(id),nonneg); 

