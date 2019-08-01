function id = mkcvindex(ns,meth,splits,strat)
% id = mkcvindex(ns,meth,splits)
% makes index vector for crosvalidation
% 
% Input:    ns - number of samples
%           meth - method leave one out ('loo'), venetian blinds ('123'),
%           continous subsets ('111'), random subsets ('rnd'); 
%           splits - number of splits/segments
% 
% Output:   id - vector with numbers from 1 to splits. 

if nargin==3; 
    strat = ones(ns,1); 
end
unstrat = unique(strat); 

switch meth
    case 'loo'
        id = 1:ns;
    case {'123','vet'}
        i = 1:splits;
        id = repmat(i',[1, ceil(ns/splits)]) ;
        id = id(1:ns);
    case '111'
        i = 1:splits;
        id = repmat(i',[ceil(ns/splits), 1]);
        id = id(1:ns);
        id = sort(id);
        id = id';
    case 'rnd'
        for j=1:length(unstrat)
            ic = ismember(strat,unstrat(j));
            nss = sum(ic); 
            i = 1:splits;
            idd = repmat(i,[ceil(nss/splits), 1]);
            idd = idd(1:nss);

            ii = randperm(nss);
            id(ic) = idd(ii);
            
        end
end

id = id';
