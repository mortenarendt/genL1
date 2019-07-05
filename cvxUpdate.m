function b = cvxUpdate(X,A,D,l)
% cd '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/MATLAB/mytools/cvx/'
% cvx_setup;
p = size(X,2);
cvx_begin
variable b(p)
minimize(norm(X - A*b') + l*norm(D*b,1))
cvx_end

% %%
% clc
% L1 = [norm(D*Bold,1) norm(D*B,1) norm(D*b,1)]
% L = [norm(X - A*Bold','fro') norm(X - A*B','fro') norm(X - A*b','fro')]
% L + 0.5*l*L1