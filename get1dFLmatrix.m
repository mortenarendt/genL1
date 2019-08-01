function FD = get1dFLmatrix(p)


% define fused lasso 
M = sparse(zeros(p));
m = -speye(p-1);
M(1:end-1,2:end) = m;
M = M + speye(p);
D = M*-1;
FD = D(1:end-1,:);
