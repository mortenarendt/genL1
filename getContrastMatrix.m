function D = getContrastMatrix(n)

M = nchoosek(1:n,2);
D = zeros(size(M,1),n);
D(sub2ind(size(D),1:size(M,1),M(:,1)')) = -1;
D(sub2ind(size(D),1:size(M,1),M(:,2)')) = 1;

