function [T v] = getASCAscores(x,d,B)

% calculate E
Xhat = d*B;
E = x - Xhat; 
[u s v] = svds(Xhat,2);
%project the E onto this model
T = x*v; 
