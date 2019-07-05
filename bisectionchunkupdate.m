function b = bisectionchunkupdate(X,a,Dorg)

[n p] = size(X);
[~, u] = genSMRupdate_dual(X,a,Dorg,Inf);
bls = genSMRupdate_dual(X,a,Dorg,0);
l2 = max(abs(u));
l2init = l2;
l1 = 0;
nbls = norm(bls);
%
c = 0;
critt = p;
% crittAct = p;
aa = 0.5;
while critt~=2 & c < maxit;
    c = c+1;
    LAMBDA = diag([ones(nFD,1)*aa; ones(p,1)*(1-aa)]);
    D = LAMBDA*Dorg;
    
    % make a new guess
    lnew = (l1+l2)/2;
    %     L(c) = lnew; lnew = l2*10
    b = colNorm(   genSMRupdate_dual(X,a,D,lnew));
    
    critt = length(unique(round(fac*(b/nbls))));
    %     crittAct = length(unique(round(fac*(b(abs(b)>eps)/nbls))))
    crittZero = sum(abs(b)<tol) / length(b);
    
    % update L
    %     if crittAct>1;
    %         l1 = lnew;
    %     elseif crittAct<1;
    %         l2 = lnew;
    %     end
    %     % update aa
    %     aa = (crittZero-0.8)*10;
    % %     aa = log(aa / (1-aa))
    %     aa = exp(aa) / (1 + exp(aa))
    %     critt = length(unique(round(fac*b)));
    if critt>2 & crittZero==0; % all over the place, more LASSO
%         l1 = lnew;
        aa = aa*0.5;
    elseif critt>2 & crittZero>0; 
        l1 = lnew;
        aa = min(aa*1.5,1)
    elseif critt<2 & crittZero==0; % flat non-zero line 
        l2 = lnew;
        aa = aa*0.5;
    elseif critt<2 & crittZero>0;  % flat zero line 
        l2= lnew;
        aa = min(aa*2.5,0.95)
    end
   
end