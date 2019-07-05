function [b u] = genL1(X,y,D,l)

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
%%
% trim off rows of D without any penalty
D = D(sum(abs(D)')>0,:);
options = optimset('Display','off');
r = qr(X,0);
% invX=r\(r'\X');
invX = sppinv(X); 
ytilde = X*invX*y; 
Dtilde = D*invX; 
H = Dtilde*Dtilde'; 
f = -Dtilde*ytilde;
m = size(D,1); 
u = quadprog(H,f,[],[],[],[],ones(m,1)*-l,ones(m,1)*l,[],options);
b = invX*(ytilde - Dtilde'*u); 


% plot(b); shg

% %%
% for ll = 1:10; 
%     l = 1/ll; 
% u = quadprog(H,f,[],[],[],[],ones(m,1)*-l,ones(m,1)*l,[],options);
% B(:,ll) = X'*a - D'*u;
% end
function invX = sppinv(X)

% invX = inv(X'*X)*X'
XtX = X'*X; 
invX = inv(XtX)*X';
% max(invX(:))

