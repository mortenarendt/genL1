%% Initialize script
clear; clc; close all;

%% Load Data
load 'RamanMeat.mat'

%% Ask about the following
% pca(auto(lgX),2);
%x = auto(lgX.data);%r = cluster(x');close all;x = x(:,r.order);
%x = [x auto(10*randn(size(x)))]; 
%clas = lgX.class{1,2};
%d = mkdesignmatrix(clas);

%% Set varibales
d = DesignTemp;
x = RamanMeat;

%% Initialize first L1manova
lambda = [0 1 0]; 
L = linspace(0,2,10);

% run first L1manova
t = cputime;
for i=1:length(L)
    tic;
    B2{i} = genL1manova(x,d,lambda/sum(lambda),L(i));
    t1 = toc;
    fprintf("Iteration: %d Completed in %f seconds \n",i,t1);
end
fprintf("Total CPU time spend on loop: %fs \n",cputime - t);

%%
for i=1:length(L)
    subplot(3,4,i);
     plot(B2{i}');
end

%%
close all; 
n = size(x,1); 
Xhatcv = x; 
for i=1:size(x,1)
    idd = true(n,1);
    idd(i) = false;
    B0 = genL1manova(x(idd,:),d(idd,:),lambda/sum(lambda),0);
    Xhatcv0(~idd,:) = d(~idd,:)*B0;
    
    B1 = genL1manova(x(idd,:),d(idd,:),lambda/sum(lambda),1.3);
    Xhatcv1(~idd,:) = d(~idd,:)*B1;
end

[~,~,v0]= svds(Xhatcv0,2); 
T0 = x*v0;

[~,~,v1]= svds(Xhatcv1,2); 
T1 = x*v1;

%% Ask about plots
%subplot(2,2,1); plotpcascores(T0,clas',[1 2],[],1); legend off; 
subplot(2,2,3); plot(v0(:,1),v0(:,2),'o'); 
% text(v0(:,1),v0(:,2),X.label{2}); 
hline(0,'k'); vline(0,'k')
%subplot(2,2,2); plotpcascores(T1,clas',[1 2],[],1); legend off; 
subplot(2,2,4); plot(v1(:,1),v1(:,2),'o'); 
% text(v1(:,1),v1(:,2),X.label{2}); 
hline(0,'k'); vline(0,'k')
