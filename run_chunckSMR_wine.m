clear; close all; clc; 
load '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/Course and teaching/Advanced Chemometrics 2013-2014/data/WineIR_Aroma.mat'

X = auto(mscorr(IR.data)); 
Y = auto(Aroma.data); 
plot(X'); shg

options = cluster('options');
options.plots = 'none';
% results = cluster(X',options);
% X = X(:,results.order); 
results = cluster(Y',options);
Y = Y(:,results.order); 

XY = corr(X,Y);
% XY = X'*Y ;
close all; 
Dx = getDfromcorrX(X,[1 1]); 
Dy = getDfromcorrX(Y,[1 1]); 
%%
k = 3; 
[A B D] = chunkSMR(XY',k,'defl',Dy,Dx); 
%%
Rhat = A*diag(D)*B';

maxcolor = max(abs(XY(:)));
figure; 
heatmap(XY',[],[],[],'Colormap', @redblue,'Colorbar', true,'MinColorValue', -maxcolor, 'MaxColorValue', maxcolor,'ShowAllTicks', true);
figure; 
heatmap(Rhat,[],[],[],'Colormap', @redblue,'Colorbar', true,'MinColorValue', -maxcolor, 'MaxColorValue', maxcolor,'ShowAllTicks', true);
%%
[A1 B1 D1] = chunkSMR(XY',k,'als'); Rhat1 = A1*diag(D1)*B1';

shg