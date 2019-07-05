function [b u fit] = genSMRupdate_dual(X,a,D,l)

% Perform generalized update of the problem: 
% |X - a*b'|fro s.t. |Db|1 < c
% via the penalty form: 
% |X - a*b'|fro + l*|Db|1
% by solving the langrange dual problem
% |X'a - D'u|2 s.t. |u|inf < l
% 
% X ~ Data (n x p)
% a ~ normalized left vector (n x 1)
% D ~ Operator matrix on b (m x p)
% l ~ L1 constraint
% 
% Morten Rasmussen 
% Marts 2016
% set parameters for the QP problem
options = optimset('Display','off');
D = sparse(D);
H = D*D';
f = -a'*X*D';
m = size(D,1); 
u = quadprog(H,f,[],[],[],[],ones(m,1)*-l,ones(m,1)*l,[],options);
b = X'*a - D'*u;

SSTOT = trace(X'*X); 
E = X - a*b';
SSE = trace(E'*E);
fit = 1 - SSE/SSTOT;

% %%
% for ll = 1:10; 
%     l = 1/ll; 
% u = quadprog(H,f,[],[],[],[],ones(m,1)*-l,ones(m,1)*l,[],options);
% B(:,ll) = X'*a - D'*u;
% end
