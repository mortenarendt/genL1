function [a d b number_breaks] = gradientFLpca(X,nc,lambda2)
% solves the Fused Lasso problem in a PCA context: 
% min_{a,b,d}: (1/2||X - a*d*b'||^2 + lambda2*|as_i - as_{i+1}|1)
% where the gradient of as (sorted version of a) is penalized. 
% This leads to grouping of the samples (a) in the scores. 

nc = 1; 
[n p] = size(X);

% some options
maxit = 30;
maxStep = 50;
thr_outer = 1e-6;
thr_inner = 1e-6; 

% initialize
[a d b] = svds(X,1);
% a = rand(n,1);
% b = rand(p,1);
c = 0;
nn = n-1;

diff = 1; 
SSEold = trace(X'*X); 
while c < maxit & diff>thr_outer
    c = c+1;
    % estimate LS
    als = X*b / sqrt(b'*b);
    als = colNorm(als); 
    [~,id] = sort(als); 
    % biSection to get the right number of clusters
    a(id) = bisectflsa(als(id),lambda2, maxStep, thr_inner);
%     [a(id), z, infor]=flsa(als(id), zeros(nn,1),0,lambda2, n,maxStep, thr_inner, 1, 6);
    a = colNorm(a);
    
    % estimate LS b
    b = X'*a / sqrt(a'*a);
    b = colNorm(b);
    d = a'*X*b;
    
    e = X - a*d*b';
    SSE = trace(e'*e);
    diff = SSEold - SSE;
    SSEold = SSE; 
end
number_breaks = length(unique(round(a*1e2)));  
disp(['Computed in ' num2str(c) ' iterations'])
disp(['SSE = ' num2str(SSE)])
% subplot(1,3,1);
% plot(SSE); shg
% subplot(1,3,2);
% plot(a(id),'.-')

% [u s v] = svds(X,1); 
% subplot(1,3,3);
% plot(t,a,'b*'); hold on; 
% plot(t,u,'r*'); hold on; 

function a = bisectflsa(als,nclus, maxStep, thr_inner)

n = length(als); 
nn = n-1; 

l1 = 0;
l2 = 1e4;

% a = flsa(als, zeros(nn,1),0,lambda2, n,maxStep, thr_inner, 1, 6);

iter = 1;
maxit = 1000;
fac = 1e2;
while iter < maxit    
%     [l1 l2]
    % make a new guess
    lnew = (l1+l2)/2;
    [a ii info] = flsa(als, zeros(nn,1),0,lnew, n,maxStep, thr_inner, 1, 6);
%     stnew = softth(t,lnew,nonneg);
    % check it
%     length(unique(a))
% unique(round(a*fac))
    if length(unique(round(a*fac)))<=nclus
        l2 = lnew;
    else
        l1= lnew;
    end
%     if length(unique(round(a*fac)))==nclus; 
    if abs(l1-l2)<thr_inner;
        disp(['#unique ' num2str(length(unique(round(a*fac))))])
        return
    end
    iter = iter+1;
    
end









