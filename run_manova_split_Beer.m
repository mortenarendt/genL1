cd '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/MATLAB/work/Sparsity_in_chemometrics/genSMR'
clear; clc; close all;
load '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/Course and teaching/ASCA belgium 2015/Data/BeerData.mat'
% pca(auto(lgX),2);
x = auto(lgX.data);%r = cluster(x');close all;x = x(:,r.order);
x = [x auto(10*randn(size(x)))]; 
clas = lgX.class{1,2};
d = mkdesignmatrix(clas);

lambda = [0 1 0]; 

L = linspace(0,2,10);
clear i B2
for i=1:length(L);
    i
    B2{i} = genL1manova(x,d,lambda/sum(lambda),L(i));    
end

for i=1:length(L);
    i
    subplot(3,4,i);
     plot(B2{i}');
     ylim([-2 2])
     
end
shg
%%
close all; 
n = size(x,1); 
Xhatcv = x; 
for i=1:size(x,1)
    idd = true(n,1);
    idd(i) = false;
    B0 = genL1manova(x(idd,:),d(idd,:),lambda/sum(lambda),0);Xhatcv0(~idd,:) = d(~idd,:)*B0;
    B1 = genL1manova(x(idd,:),d(idd,:),lambda/sum(lambda),1.3);Xhatcv1(~idd,:) = d(~idd,:)*B1;
end

[u s v0]= svds(Xhatcv0,2); T0 = x*v0;
[u s v1]= svds(Xhatcv1,2); T1 = x*v1;

% B0 = genL1manova(x,d,lambda/sum(lambda),0);    
% [T0 v0] = getASCAscores(x,d,B0);
% B1 = genL1manova(x,d,lambda/sum(lambda),1.3);    
% [T1 v1] = getASCAscores(x,d,B1);

%%
subplot(2,2,1); plotpcascores(T0,clas',[1 2],[],1); legend off; 
subplot(2,2,3); plot(v0(:,1),v0(:,2),'o'); 
% text(v0(:,1),v0(:,2),X.label{2}); 
hline(0,'k'); vline(0,'k')
subplot(2,2,2); plotpcascores(T1,clas',[1 2],[],1); legend off; 
subplot(2,2,4); plot(v1(:,1),v1(:,2),'o'); 
% text(v1(:,1),v1(:,2),X.label{2}); 
hline(0,'k'); vline(0,'k')



shg








