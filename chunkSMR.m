function [A B D] = chunkSMR(X,k,method,D1,D2)

if nargin==2;
    method = 'als';
end

if nargin<4;
    D1 = []; D2 = [];
end

switch method
    case 'als'
        diff = 1;
        conv = 1e-5;
        maxit = 10;
        c = 0;
        
        % initialize
                 [A D B] = svds(X,k);D = diag(D);
%         A = zeros(size(X,1),k); B = zeros(size(X,2),k); D = zeros(k,1);
        
        SSEold = trace(X'*X);
        while diff>conv & c<maxit
            % all components simultanously
            for i=1:k;
                idout = true(k,1);
                idout(i) = false;
                Xtilde = X - A(:,idout)*diag(D(idout))*B(:,idout)';
                [A(:,i) B(:,i) D(i)] = chunkSMR_onecomp(Xtilde,D1,D2,A(:,i));
            end
            
            E = X  - A*diag(D)*B';
            SSE = trace(E'*E);
            %         sse(c) = SSE;
            diff = SSEold - SSE;
            SSEold = SSE;
        end
        
    case 'defl'
        for i=1:k;
            [a b d] = chunkSMR_onecomp(X,D1,D2);
            A(:,i) = a;
            B(:,i) = b;
            D(i) = d;
            X = X - a*d*b';
        end
end

function [a b d] = chunkSMR_onecomp(X,D1,D2,a)
SStot = trace(X'*X);
SSEold = SStot;

if nargin==3;
    % initialize
    [a s] = svds(X,1);
end

diff = 1;
conv = 1e-5;
maxit = 100;
c = 0;

while diff>conv & c<maxit
    c = c+1;
    % update b
    b = colNorm(chunkupdate(X,a,D2));
    
    % update a
    a = colNorm(chunkupdate(X',b,D1));
    d = a'*X*b;
    E = X  - a*d*b';
    SSE = trace(E'*E);
    sse(c) = SSE;
    diff = SSEold - SSE;
    SSEold = SSE;
end
