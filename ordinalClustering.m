function [a b d SSE] = ordinalClustering(x,k,l,nonneg)
%%%%%% input 
% x     - data (n,p)
% k     - # components (1,1)
% l     - penalty (1,1)

% x = auto(X.data);
% l = 3;

% some settings
conv = 1e-6;
maxit = 100; 

[n p] = size(x);
D = get1dFLmatrix(n);
% initialize
[u d v] = svds(x,k); d = diag(d); 
sng = sign(sum(corr(x',v)));
a = u*diag(sng); b = v*diag(sng);

% a = colNorm(rand(n,k)); b = colNorm(rand(p,k)); d = ones(1,k); 

SSEold = trace(x'*x);
c = 0; 
diff = 1;
while (diff>conv & c<maxit) %| c<3
    c = c+1;
    % all components simultanously
    for i=1:k;
        idout = true(k,1);
        idout(i) = false;
        Xtilde = x - a(:,idout)*diag(d(idout))*b(:,idout)';
        [a(:,i), d(i), b(:,i)] = ordinalcomponentupdate(Xtilde,a(:,i),b(:,i),D,l,nonneg);
    end
    
    E = x  - a*diag(d)*b';
    SSE = trace(E'*E);
    diff = SSEold - SSE;
    SSEold = SSE;
end

