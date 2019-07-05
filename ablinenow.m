function [r p h] = ablinenow(x,y)

if nargin==0; 
x = get(get(gca,'children'),'XData');
y = get(get(gca,'children'),'YData');
end

if size(x,1)==1; 
    x = x'; 
end

if size(y,1)==1; 
    y = y'; 
end


id1 = sum(isnan([x y]'))==0; 
id2 = sum(isinf([x y]'))==0; 
id = id1 & id2; 

x = x(id); 
y = y(id); 
n = length(x);

X = [ones(n,1) x]; 
b = pinv(X)*y; 
hold on; 
h = abline(b(2),b(1)); 
hold off

[r p ] = nancorr(x,y);



