function l=searchLuc(t,sumabst,nonneg)

% seach for the lambda penalty for fixing the constraint: 
% sum(abs(t)) < sumabst via soft-thressholding

if norm(t,2)==0 || sum(abs(t))<=sumabst
    l=0;
    return
end
l1 = 0;
l2 = max(abs(t))-1e-4;

iter = 1;
maxit = 1000;
while iter < maxit
    % make a new guess
    lnew = (l1+l2)/2;
    stnew = softth(t,lnew,nonneg);
    % check it
    if sum(abs(stnew))<sumabst
        l2 = lnew;
    else
        l1= lnew;
    end
    if (l2-l1)<1e-5
        l=lnew;
        return
    end
    iter = iter+1;
end

