function b = chunkupdate(X,a,Dorg)

lambda1 = 1; % fused lasso penalty
lambda2 = 1; % lasso penalty
maxit = 10;
fac = 1e2; % rounding factor of similarity.
[n p] = size(X);

if nargin==2 | isempty(Dorg);
    % define fused lasso and lasso penalty
    M = zeros(p);
    m = -eye(p-1);
    M(1:end-1,2:end) = m;
    M = M + eye(p);
    D = M*-1;
    FD = D(1:end-1,:);
    Dorg = [FD; eye(p)];
end
nFD = p-1; 
LAMBDA = sparse(diag([ones(nFD,1)*lambda1*0.5; ones(p,1)*lambda2*0.5])); 
Dorg = sparse(Dorg); 
D = LAMBDA*Dorg;

b = gridchunkupdate(X,a,Dorg,10);
% b = bisectionchunkupdate(X,a,Dorg);

% [~, u] = genSMRupdate_dual(X,a,D,Inf);
% bls = genSMRupdate_dual(X,a,D,0);
% nbls = norm(bls);
% 
% aa = linspace(0,1,10);
% LL = linspace(0,max(abs(u)),10);
% BB = []; b = []; 
% for i=1:length(aa); 
%     LAMBDA = diag([ones(nFD,1)*aa(i); ones(p,1)*(1-aa(i))]);
%     D = LAMBDA*Dorg;
%     for ii=1:length(LL); 
%         bb = genSMRupdate_dual(X,a,D,LL(ii));
%         critt = length(unique(round(fac*(bb/nbls))));
%         b(ii,:) = [bb; aa(i); LL(ii); critt]';
%         
%     end
%     BB = [BB; b];     
% end
% %
% BBx = dataset(BB(:,1:end-3)); 
% BBx.axisscale{1,1} = BB(:,end-1); BBx.axisscalename{1,1} = 'L';
% BBx.axisscale{1,2} = BB(:,end-2); BBx.axisscalename{1,2} = 'a';
% BBx.axisscale{1,3} = BB(:,end); BBx.axisscalename{1,3} = 'Critt';
% BBx.class{1,3} = BB(:,end); BBx.classname{1,3} = 'Critt';
% pca(auto(BBx),2)
% %%
% close all; 
% plot(BBx.axisscale{1,1},BBx.axisscale{1,2},'o'); xlabel('L'); ylabel('a');
% hold on; 
% idd = BBx.axisscale{1,3}==2; 
% plot(BBx.axisscale{1,1}(idd),BBx.axisscale{1,2}(idd),'*', 'markersize',9);
% hold off; 
% text(BBx.axisscale{1,1},BBx.axisscale{1,2},num2str(BBx.axisscale{1,3}')); 
% ylim([-.1 1.1]); xlim([-.1 max(LL)+0.1]); shg
% figure; 
% plot(BBx.data(idd,:)'); 
% %%
% 
% get maximum lambda (max(abs(u) unconstrained)


% % update LS b with these two 
% nonzero = abs(b)>eps;
% b = zeros(length(b),1); 
% b(nonzero) = mean(X(:,nonzero)'*a); 
