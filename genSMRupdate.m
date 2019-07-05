function b = genSMRupdate(X,a,D,l)

% Perform generalized update of the problem: 
% |X - a*b'|fro s.t. |Db|1 < l
% 
% X ~ Data (n x p)
% a ~ normalized left vector (n x 1)
% D ~ Operator matrix on b (p x p)
% l ~ L1 constraint
%
% Method: 
% 1) transform the problem in to a normal LASSO problem, by setting: z = Db
% (D needs to be intertible for this to work)
% 2) Solve by softthresholding of LeastSquares estimates
% 3) Backtransformation

% Calculate the inverse of D 
Dinv = pinv(D);

% LS update of z
z = D*X'*a / (a'*a); 

% softthresshold z by relevant lambda
%%% find relevant penalty
ll = searchLuc(z,l,0);
%%% softthresshold by this one
z = softth(z,ll);

% backtransformation 
b = Dinv*z;
