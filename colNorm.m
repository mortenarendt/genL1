function A = colNorm(A)

A = A*diag(1./sqrt(diag(A'*A)));
A(isnan(A)) = 0; 
