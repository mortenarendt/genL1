cd '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/MATLAB/work/Sparsity_in_chemometrics/genSMR'
clear; clc; close all;
load '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/Course and teaching/ASCA belgium 2015/Data/BeerData.mat'
% pca(auto(lgX),2);
% x = auto(lgX.data);r = cluster(x');close all;x = x(:,r.order);d = mkdesignmatrix(lgX.class{1,2});
load '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/Course and teaching/ASCA belgium 2015/Data/MeatRaman.mat'
[~,id1] = min(abs(Xpp.axisscale{2} - 700))
[~,id2] = min(abs(Xpp.axisscale{2} - 1100))

% select temperatures
tempsel = [52 54 62 70];
tempid = ismember(Xpp.class{1,1},tempsel)
timeid = ismember(Xpp.class{1,2},[2])
% pca(mncn(Xpp),2)
% X = Xpp(tempid & timeid,id1:id2); 
X = Xpp(tempid & timeid,:); 

% pca(mncn(x),2)
%
% x = mncn(Xpp.data(:,1000:1:1500));
plot(X'); shg
% plot(x.data'); shg
% pca(x,2); shg
d = mkdesignmatrix(X.class{1});
x = mncn(X.data); 
%%
lambda = [1 1 1]; 
%%%%%%%%%%%%%%

L = linspace(0,300,10);
clear i B2
for i=1:length(L);
    tic
    i
    B2{i} = genL1manova(x,d,lambda/sum(lambda),L(i));    
    toc
end

for i=1:length(L);
    i
    subplot(3,4,i);
     plot(B2{i}');
%     ylim([-2 2])
end

save penalizedManovaMeatRaman.mat
%%
clear; 
load penalizedManovaMeatRaman.mat
close all; 
pt = 30; 
legFS = 14; 
% select a single solution 
xax = X.axisscale{2};
leg = [repmat('T = ',4,1) num2str(tempsel') repmat('^oC',4,1)];

figure
B = B2{9}; 
id = unique([1:pt:size(B,2) find(sum(abs(B))>0.1)])
plot(xax(id),30*B(:,id)'); hold on; 
plot(xax(id),X.data(1:7:end,id) - min(X.data(:)),'color',ones(1,3)*0.5); shg
h = legend(leg,'location','SouthWest')
set(h, 'FontSize',legFS); 
set(gca,'xdir','reverse');
axis tight; 
hold off;
xlabel('Raman shift (cm^{-1})','interpret','tex')
set(gca,'ytick',0)
% matlab2tikz('paper/Fig4B.tex','width','3in','height','3in'); close all; 

print('paper/Fig4B_v2','-depsc'); close all; 
% close all; 
figure
B = B2{1}; 
plot(xax(1:pt:end),30*B(:,1:pt:end)'); hold on; 
plot(xax(1:pt:end),X.data(1:2:end,1:pt:end) - min(X.data(:)),'color',ones(1,3)*0.5); shg
h2 = legend(leg,'location','SouthWest')
set(h2,'FontSize',legFS)
set(gca,'xdir','reverse');
axis tight; 
hold off;
xlabel('Raman shift (cm^{-1})','interpret','tex')
set(gca,'ytick',0)
% matlab2tikz('paper/Fig4A.tex','width','3in','height','3in'); close all; 
print('paper/Fig4A_v2','-depsc'); close all; 




