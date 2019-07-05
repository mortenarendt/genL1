function P = simulateChunkloading(p,k,nactive)


for i=1:k
    PX = zeros(p,1);
    
    cnt = round(rand(1)*(p-nactive)) + nactive/2;
    st = max(round(cnt - nactive/2),1);
    sl = min(round(cnt + nactive/2),p);
    
    PX(st:sl) = 1;
    P(:,i) = PX;
end

