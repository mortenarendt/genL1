function [B,b,u] = genL1manova(x,d,lambda,L)

% INPUT: 
%     x~(n,p)     - Data
%     d~(n,m)     - design
%     lambda      - "simplex" [lasso penalty, group contrast penatly, spectral FL penalty, ]
%     L           - overall penalty
% OUTPUT: 
%     B~(m,p)     - Regression Coefficients (in original shape)
%     b~(mp,1)    - Regression Coefficients vectorized
%     u~(size(D,2),1) - dual solution

% get dimentions
[~,p] = size(x);
m = size(d,2);

% define the penalty matrices
% FL on the feature dimention
e = ones(p,1);
Dfl1 = spdiags([-e  e], [0,1], p-1, p);

% Contrasts on design
Dfl2  = getContrastMatrix(m);

% vectorize X and design
vx = x(:);
vd = kron(speye(p),sparse(d));

% set penalty matrix accordingly
vDfl2 = kron(speye(p),Dfl2);
vDfl1 = kron(Dfl1,speye(m));
D = [vDfl1*lambda(3); vDfl2*lambda(2); speye(p*m)*lambda(1)];

% Calculating result
[b,u] = genL1(vd,vx,D,L);
B = reshape(b,m,p);
