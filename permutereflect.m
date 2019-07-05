function res = permutereflect(T,R,ref)

ncomp = size(T,2);

% number of permutations
% permutation matrix (P)
P = perms(1:ncomp);

np = size(P,1);

for i=1:np;
    cc(i) = tuckersCC(T,R(:,P(i,:)),ref);
end

res.cc = cc;
res.P = P;
I = eye(ncomp);
[m id] = max(cc);
p = P(id,:);
Q = I(p,:);

% permute optimaly
RR  = R(:,p);

clear cc;
if ref==1;
    % reflect
    for i=1:ncomp;
        cc(i) = corr(T(:,i),RR(:,i));
    end
    S = ones(ncomp,1)*sign(cc);
    res.Q = Q.*S';
else
    res.Q = Q;
end

function c = tuckersCC(T,Tb,ref)
p = size(T,2);
for i = 1:p;
    t = T(:,i);
    nt = sqrt(t'*t);
    tb = Tb(:,i);
    ntb = sqrt(tb'*tb);
    c(i) = t'*tb / (nt*ntb);
end
if ref==1;
    c = sum(abs(c));
else
    c = sum(c);
end

