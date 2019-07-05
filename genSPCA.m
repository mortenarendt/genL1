function results = genSPCA(X,D,k,l,init,method)

if ismember(method,'cvx')
    try
    cvx_setup
    catch
        disp('Download and install the cvx lib')
        return
    end
end

[n p] = size(X);

if ismatrix(D) % only constraint on second mode
    DD{1} = eye(n);
    DD{2} = D;
    D = DD;
end

if length(l)==1; % only constraint on second mode
%     L = [Inf l];
L = [0 l];
    l = L;
end

results.lambda = l;
SSETOT = norm(X,'fro')^2;

if ischar(init)
    switch init
        case 'svd'
            % Initialize by SVD
            [A,s,B] = svds(X,k);
        case 'rnd'
            % initialize by random numbers
            A = colNorm(randn(n,k));
            B = colNorm(randn(p,k));
    end
else
    A = init{1};
    B = init{2};
end

for i=1:size(B,2)
    z = D{2}*B(:,i);
    ll = searchLuc(z,l(2),0);
    %%% softthresshold by this one
    B(:,i)= pinv(D{2})*softth(z,ll);
end

results.init{1} = A;
results.init{2} = B;
results.sstot = SSETOT;
Fit = 0;
conv =1e-6;
diff = 1;
c = 0;
compUD = 1:k;
while diff > conv
    c = c+1;
%     % Update A
%     [A cont(1,:)] = genSMRupdateM(X',B,A,D{1},l(1),method);
%     % normalize A
%     A = colNorm(A); 
%     % Update B (constrained)
%     [B cont(2,:)] = genSMRupdateM(X,A,B,D{2},l(2),method);
%     
    for j=1:k
        i = compUD(j);
        noti = true(1,k);
        noti(i) = false;
        Xtilde = X - A(:,noti)*B(:,noti)';
        % Update A
        [A(:,~noti) cont(1,~noti)] = genSMRupdateM(Xtilde',B(:,~noti),A(:,~noti),D{1},l(1),method);
        % normalize A
        A = colNorm(A); 
        % Update B (constrained)
        [B(:,~noti) cont(2,~noti)] = genSMRupdateM(Xtilde,A(:,~noti),B(:,~noti),D{2},l(2),method);
    end
    
    Fitold = Fit;
    SSE = norm(X - A*B','fro')^2;
    Fit = (SSETOT - SSE)/SSETOT;
    diff = (Fit - Fitold);
    FIT(c) = Fit;
    sse(c) = SSE;
end

expvar = calexpvar(X,A,B);
results.Fit = FIT;
results.sse = sse;
results.expvar = expvar;
results.loads{1} = A;
results.loads{2} = B;
results.constraint_active = cont;


%%%%%%%%% INTERNAL FUNCTIONS
function expvar = calexpvar(data,scores,loads)

idnan = isnan(data); 
SStot = data(idnan==0);
SStot = SStot'*SStot;
modelX = scores*loads';
SSmod = modelX(idnan==0);
SSmod = SSmod'*SSmod;
EVallpc = SSmod/SStot;

ncomp = size(scores,2);

[ns,nv] = size(data); 

for i=1:ncomp; 
    modelX = scores(:,i)*loads(:,i)';
    SSmod = modelX(idnan==0);
    SSmod = SSmod'*SSmod;
    EV = SSmod/SStot;
%     exp(i) = EV;
    exppc(i) = EV;
    
    % incremental
    modelX = scores(:,1:i)*loads(:,1:i)';
    SSmod = modelX(idnan==0);
    SSmod = SSmod'*SSmod;
    EV = SSmod/SStot;
%     exp(i) = EV;
    exppc_incremental(i) = EV;
    
    idd = false(ncomp,1); 
    idd(i) = true; 
     
    modelX = scores(:,~idd)*loads(:,~idd)';
    SSmod = modelX(idnan==0);
    SSmod = SSmod'*SSmod;
    EV = SSmod/SStot;
    margpc(i) = EVallpc - EV; 
    for ii=1:nv; 
        modelXvar = modelX(:,ii); 
        Xvar = data(:,ii); 
        resXvar = Xvar - modelXvar; 
        
        SSe = resXvar(idnan(:,ii)==0); 
        SSe = SSe'*SSe; 
        
%         SSmodvar = modelXvar(idnan(:,ii)==0);
%         m = SSmodvar; 
%         SSmodvar = SSmodvar'*SSmodvar;
        SStotvar = Xvar(idnan(:,ii)==0); 
%         t = SStotvar;
        SStotvar = SStotvar'*SStotvar;
        
        
        EVvar(i,ii) = 1 - SSe/SStotvar;
    end
end

EVvar(EVvar==Inf) =0;
% expvar.total = exp; 
expvar.marginal = margpc; 
expvar.variable = EVvar; 
expvar.prpc = exppc;
expvar.incremental = exppc_incremental; 

