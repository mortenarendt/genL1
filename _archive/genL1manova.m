function [B b u] = genL1manova(x,d,lambda,L)

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
[n p] = size(x);
m = size(d,2);

% define the penalty matrices
% FL on the feature dimention
Dfl1 = get1dFLmatrix(p);
% Contrasts on design
Dfl2  = getContrastMatrix(m);

% vectorize X and design
vx = vec(x);
% vd = kron(eye(p),d);
vd = kron(sparse(eye(p)),sparse(d));

% set penalty matrix accordingly
vDfl2 = kron(eye(p),Dfl2);
vDfl1 = kron(Dfl1,eye(m));

D = [vDfl1*lambda(3); vDfl2*lambda(2); eye(p*m)*lambda(1)];

% remove zero rows from D 
D = D(sum(abs(D'))>0,:);

D = sparse(D); 
vd = sparse(vd); 

[b u] = genL1(vd,vx,D,L);
% % % LS update of active b's 
% % isZero = abs(b)>1e-10;
% % b = zeros(length(b),1); 
% % b(isZero) = sppinv(vd(:,isZero))*vx;
B = reshape(b,m,p);


function invX = sppinv(X)

% invX = inv(X'*X)*X'
XtX = X'*X; 
invX = inv(XtX)*X';
% max(invX(:))