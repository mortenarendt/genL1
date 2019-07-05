function b = gridchunkupdate(X,a,Dorg,npoints)
%%
[n p] = size(X); 
nFD = p-1; 
fac = 1e2; % rounding factor of similarity.
tol = 1e-4; 
% set max on the L scale
[~, u] = genSMRupdate_dual(X,a,Dorg,Inf);
bls = colNorm(genSMRupdate_dual(X,a,Dorg,0));
nbls = norm(bls);

% construct grid
aa = linspace(0,1,npoints);
LL = linspace(0,max(abs(u)),npoints);
BB = []; bb = [];
for i=1:length(aa);
    LAMBDA = diag([ones(nFD,1)*aa(i); ones(p,1)*(1-aa(i))]);
    D = LAMBDA*Dorg;
    for ii=1:length(LL);
        [b u fit] = genSMRupdate_dual(X,a,D,LL(ii));
        b = colNorm(b);
        
        % evaluate in terms of 1) number of zeros and 2) fusion in the
        % actives
        crittAct = length(unique(round(fac*(b(abs(b)>tol)/nbls))));
        crittZero = sum(abs(b)<tol) / length(b);
        
        % collect results.
        bb(ii,:) = [b; aa(i); LL(ii); crittAct; crittZero; fit]';
        
    end
    BB = [BB; bb];
end
%%
% find optimal solution. 
BB(:,1:p) = colNorm(BB(:,1:p)')';


idd1 = BB(:,end-2); 
idd2 = BB(:,end-1);
idd3 = BB(:,end);

% obj = idd1/range(idd1) - idd2/range(idd2)*0.1;
obj = -idd3;
obj(idd1==0 | idd1 > 2 | idd2<0.2 ) = NaN; 
% plot(obj,'*'); shg
%
[~,idd] = nanmin(obj);
b = BB(idd,1:p)'; 
% plot(BB(idd,1:p))
% %%
% idd1==1;
% % idd2(idd1~=1) = 0; 
% [~,idd] = max(idd2)
% b = BB(:,1:p);
% close all; 
plot(b'); shg
% plot(idd1,idd2,'o'); shg
% % plot(BB(idd1==1,1:179)'); shg
% % update LS b with these two
% nonzero = abs(b)>eps;
% b = zeros(length(b),1);
% b(nonzero) = mean(X(:,nonzero)'*a);
