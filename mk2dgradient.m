function D = mk2dgradient(nx,ny)
% Make D
% nx = length(x);
% ny = length(y);
% nx = 5;
% ny = 6;
M = zeros(nx);
m = -eye(nx-1);
M(1:end-1,2:end) = m;
M = M + eye(nx);
M = M(1:end-1,:)';
d = M'*-1;
d0 = zeros(size(d));
D = [];
for i=1:ny;
    D = [D; [repmat(d0,1,i-1) d repmat(d0,1,ny-i)]];
end

d0 = zeros(ny);
d1 = eye(ny);
for i = 1:nx-1;
    
    D = [D; [repmat(d0,1,i-1) d1 -d1 repmat(d0,1,nx-i-1)]];
end
