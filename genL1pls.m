function [T P W B] = genL1pls(x,y,lv,D,options)

% options   - vector: option(1): LS update on the active elements of w Yes
% = 1 / No = 0 (default = 0)

yorg = y;
[m nx] = size(x);
ny = size(y,2);
thr = 1e-7;
%
for i = 1:lv
    % calculate genL1 penalized weights
    [rr u] = genSMRupdate_dual(x,y,D,1);
    % check if all == 0
    icact = abs(rr)>thr;
    if sum(icact)==0
        disp('No active variables')
    else
        rr = colNorm(rr);
        if options(1)==1
        % update rr with only the active features
        rr(icact) = x(:,icact)'*y;
        rr = colNorm(rr);
        end
    end
    qq = 1;
    
    %calc loads and scores for this component
    tt = x*rr;
    normtt = norm(tt);
    NT(i) = normtt;
    %   tt = tt/normtt;
    %   rr = rr/normtt;
    pp = (tt'*x)';
    
    
    qq = y'*tt;
    uu = y*qq;
    vv = pp;
    %     if i > 1
    %         vv = vv - basis*(basis'*pp);
    %         uu = uu - xscrs*(uu'*xscrs)';
    %     end
    vv = vv/norm(vv);
    
    % deflate y
    y = (eye(m) - tt*pinv(tt))*y;
    
    % deflate X
%     x = (eye(m) - tt*pinv(tt))*x;
    
    wts(:,i)   = rr;           %r  x-block weights
    xscrs(:,i) = tt;           %t  x-block scores
    xlds(:,i)  = pp;           %p  x-block loadings
    ylds(:,i)  = qq;           %q  y-block loadings
    yscrs(:,i) = uu;           %u  y-block scores
    basis(:,i) = vv;           %v  basis of x-loadings??
end

T = x*wts;
P = xlds;
W = wts;
for i = 1:lv;
    B(:,i)  = W(:,1:i)*pinv(T(:,1:i))*yorg;
end
