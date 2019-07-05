clear; clc;
cd '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/MATLAB/work/Sparsity_in_chemometrics/genSMR'
load HPLCforweb.mat
% remove sample 119
X = HPLCforweb([1:118 120],:);
% do baseline correction
X = wlsbaseline(X);
Xorg = X;
X.data = mncn(X.data);
X.data = X.data / norm(X.data);

[~,id] = sort(X.class{1,1});
X = X(id,:);

L = linspace(0,2,20); 
for i=1:length(L); 
% make ordinal model
[a b d sse] = ordinalClustering(X.data,2,L(i),0);
fac = 1000; 
una(i) = length(unique(round(a(:,1)*fac))/fac); 
end
%%
[a b d sse] = ordinalClustering(X.data,2,1.1,0);

sng = [1 -1];
a = a*diag(sng); b = b*diag(sng);
% make unconstrained model
[a2 d2 b2] = svds(X.data,2);

sng = sign(diag(a'*a2));
a2 = a2*diag(sng);
b2 = b2*diag(sng);

%% plot it
close all;
% loadingplot
leg = {'Unconstrained';'Ordinal factor model'};
i = 1;
rt = linspace(0,600,size(X,2));
plot([b2(:,i) b(:,i)]); hold on;
% plot(0.005*Xorg.data','color',ones(1,3)*0.3)
axis tight;
legend(leg,'location','southeast')
xlabel('Retention Time')

height = 10; width = 10; setfiguresize
% matlab2tikz('paper/Fig22A.tex'); close all; 
% print('paper/Fig22A','-depsc'); close all;

%%
close all;
xax = 0;
% score plot
i = 1;
A = [a2(:,i) a(:,i)];
clss = X.class{1,1};
leg2 =X.classlookup{1,1}(2:end,2);
lg = repmat(leg2',2,1); lg = [char(lg(:)) repmat(' ',6,1) char(repmat(leg,3,1))];
col = colormap(parula(4));
col = col([1:2 3],:);
for ii = 1:3;
    id= clss==ii;
    mi = max(xax)+2;
    xax = mi:(mi+sum(id)-1);
    
    plot(xax,A(id,1),'o','color',col(ii,:),'markerfacecolor',col(ii,:)); hold on;
    plot(xax,A(id,2),'*','color',col(ii,:),'markerfacecolor',col(ii,:)); hold on;
end
fac = 1000; 
una = unique(round(a(:,i)*fac))/fac; 
length(una)
hline(una,':k')
axisalmosttight
xlim([-1 max(xax)+3])
hold off;
legend(lg)
set(gca,'xtick',[]); shg
xlabel('Sample index')
ylabel('Scores component 1')
height = 10; width = 20; setfiguresize
% matlab2tikz('paper/Fig22B.tex'); close all; 
% print('paper/Fig22B','-depsc'); close all;

%%

close all;
xax = 0;
% score plot
i = 1;
A = [a2(:,i) a(:,i)];
clss = X.class{1,1};
leg2 =X.classlookup{1,1}(2:end,2);
lg = repmat(leg2',2,1); lg = [char(lg(:)) repmat(' ',6,1) char(repmat(leg,3,1))];
col = colormap(parula(4));
col = col([1:2 3],:);
[~,idd] = sort(A(:,1)); 
clss = clss(idd); 
A = A(idd,:); 
for ii = 1:3;
    id= clss==ii;
    mi = max(xax)+2;
    xax = mi:(mi+sum(id)-1);
    
%     plot(xax,A(id,1),'o','color',col(ii,:),'markerfacecolor',col(ii,:)); hold on;
    plot(A(id,1),A(id,2),'-*','color',col(ii,:),'markerfacecolor',col(ii,:)); hold on;
end
hold off; 
plot(A(:,1),A(:,2),'-*','color','r'); 

fs = 15;
xlabel('Least Squares Estimates','fontsize',fs)
ylabel('Ordinal Penalized Estimates','fontsize',fs)


