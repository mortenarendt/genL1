function [B constraintactive] = genSMRupdateM(X,A,B,D,l,method)

% Perform generalized update of the problem:
% |X - A*B'|fro s.t. |DBi|1 < l (where Bi er the individual vectors of B)
%
% X ~ Data (n x p)
% A ~ left factor matrix (n x k)
% B ~ right factor matrix (p x k)
% D ~ Operator matrix on b (p x p)
% l ~ L1 constraint (1 x 1) (same for all vectors)
%
% Method:
% 1) transform the problem in to a normal LASSO problem, by setting: z = DB
% (D needs to be intertible for this to work)
% 2) Solve by softthresholding of LeastSquares estimates
% 3) Backtransformation

%%% initially check whether the LS solution are within the contrained set.
Bold = B;
% init fit
SSEinit = norm(X - A*B','fro');

Bls = X'*A*pinv(A'*A);
k = size(Bls,2);
for i=1:k;
    CS(i) = norm(D*Bls(:,i),1);
end

constraintactive = zeros(1,k);
    % if sum(CS>l)==0
if l==0
    B = Bls;
    return
end

% Calculate the inverse of D
Dinv = pinv(D);
k = size(A,2);

% find components violating constraints
% compUD = find(CS>l);
compUD = 1:k; 
% update over the components violating the constraint
for j=1:length(compUD);
    i = compUD(j);
    noti = true(1,k);
    noti(i) = false;
    
    Xtilde = X - A(:,noti)*B(:,noti)';
    
    switch method
        case 'lasso'
            % LS update of z
            z = D*Xtilde'*A(:,i) / (A(:,i)'*A(:,i));
            
            % softthresshold z by relevant lambda
            %%% find relevant penalty
            %     ll = searchLuc(z,l,0);
            ll = l;
            constraintactive(i) = (ll>eps)+0;
            %%% softthresshold by this one
            z = softth(z,l);
            
            % backtransformation
            B(:,i) = Dinv*z;
        case 'cvx'
            B(:,i) = cvxUpdate(Xtilde,A(:,i),D,2*l);
        case 'dual'
            B(:,i) = genSMRupdate_dual(Xtilde,A(:,i),D,2*l);
        case 'dualpath'
            B(:,i) = genSMRupdate_dualpath(Xtilde,A(:,i),D,2*l);
    end
    % check fit
    SSE = norm(X - A*B','fro');
%     if (SSE - SSEinit)>0
%         B(:,i) = Bold(:,i);
% %         %%
% %         a = linspace(0,1,30);
% %         for jj=1:length(a)
% %             BB = (1-a(jj))*Bold + a(jj)*B;
% %             nBB(jj) = norm(D*BB,1);
% %             SSEb(jj) = norm(X - A*BB','fro');
% %         end
% %         subplot(1,2,1); plot(a,nBB,'o-'); shg
% %         subplot(1,2,2); plot(a,SSEb,'o-'); shg
%         
%         %%
%     else
%         Bold(:,i) = B(:,i);
%     end
end
