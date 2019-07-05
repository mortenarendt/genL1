function D = getDfromcorrX(X,lambda)

pX = size(X,2); 
% cX = nancorr(X);
cX = corr(X); 
m = diag(cX(1:end-1,2:end));
m(m<0) = 0; 
length(m)
% m = m - min(m); 
% m = m/max(m); 
% m = m + 0.1; 
% scale to max 1
m = m / mean(m);
% m = sqrt(m); 
% m = m.^2;
% m = m / mean(m);
M = zeros(pX-1,pX);
M(:,2:end) = diag(m);
M(:,1:end-1) = M(:,1:end-1) + diag(-m);
D = [M*lambda(1); eye(pX)*lambda(2)];
