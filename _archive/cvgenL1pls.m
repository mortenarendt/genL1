function results = cvgenL1pls(x,y,lv,lambda,npt,opt)

if numel(npt)==1; 
    npt = [npt npt]; 
end
% get penalty matrix
[n p]= size(x);
flD = get1dFLmatrix(p);
lD = eye(p);
thr = 1e-7;
% set range for lamdba
% LL(1,:) = exp(linspace(log(thr),log(lambda(1)),npt));
% LL(2,:) = exp(linspace(log(thr),log(lambda(2)),npt));
LL{1} = linspace(thr,lambda(1),npt(1));
LL{2}= linspace(thr,lambda(2),npt(2));
% LL = fliplr(1./LL);
% LL = LL(:,2:end); 
% generate cvindex
nsplit = round(sqrt(n));
cvID = mkcvindex(n,'rnd',nsplit);
YHAT = []; design = []; BB = []; 
for i = 1:npt(1)
    i
    for ii=1:npt(2)
        D = [flD*LL{2}(ii); lD*LL{1}(i)];
        % CV loop
        yhat = NaN(n,lv);
        for j = 1:nsplit
            ic = cvID~=j;
            % run model
            [T P W B] = genL1pls(x(ic,:),y(ic),lv,D,opt);
            % apply model
            yhat(~ic,:) = x(~ic,:)*B;
        end
        YHAT  = [YHAT yhat]; 
        design = [design [1:lv; repmat([LL{1}(i) LL{2}(ii)]',1,lv)]];
        
        [T P W B] = genL1pls(x,y,lv,D,opt);
        BB = [BB B];
    end
end

designlb = {'#lv';'LASSO penalty'; 'FusedLASSO penalty'}; 

E = YHAT - y*ones(1,size(YHAT,2));
RMSE = sqrt(diag(E'*E) / n); 

results.RMSE = RMSE;
results.Yhat = YHAT; 
results.B = BB; 
results.design = design;
results.designlb = designlb;
results.cvid = cvID; 
results.y = y; 
results.x = x; 
