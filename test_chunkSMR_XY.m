cd '/Users/mortenarendtrasmussen/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/MATLAB/work/Sparsity_in_chemometrics/genSMR'
clear; clc; close all;
rng(1565)
n1 = 100; 
n = 5000;
pX = 40;
pY = 50;
lambda1 = 1; 
lambda2 = 1; 
k = 3;
s = 3;
% on common sparse component + noise
PX = simulateChunkloading(pX,k,round(pX*0.2))*diag(sign(rand(k,1)-0.5));
PY = simulateChunkloading(pY,k,round(pY*0.2))*diag(sign(rand(k,1)-0.5));
T = randn(n,k);
Sx = T*PX'; 
Ex = randn(n,pX)*s;
Sy = T*PY'; Ey = randn(n,pY)*s;
X = Sx + Ex;
Y = Sy + Ey;
trace(Sx'*Sx) / trace(Ex'*Ex)
trace(Sy'*Sy) / trace(Ey'*Ey) 

%
% R0 = corr(T*PX'.T*PY'); 
Xt = X((n1+1):end,:); X = X(1:n1,:);
Yt = Y((n1+1):end,:); Y = Y(1:n1,:);
%
% options = cluster('options');
% options.plots = 'none';
% results = cluster(X',options);
% X = X(:,results.order); 
% Xt = Xt(:,results.order); 
% PX = PX(results.order,:); 
% results = cluster(Y',options);
% Y = Y(:,results.order); 
% Yt = Yt(:,results.order); 
% PY = PY(results.order,:); 
%
Roracle = corr(T*PX',T*PY');
Roracle(isnan(Roracle)) = 0; 
R = corr(X,Y);
Rt = corr(Xt,Yt);

% make D1 and D2 from structure in X and Y
% Dx = getDfromcorrX(R',[lambda1 lambda2]); 
% Dy = getDfromcorrX(R,[lambda1 lambda2]); 
Dx = getDfromcorrX(X,[lambda1 lambda2]); 
Dy = getDfromcorrX(Y,[lambda1 lambda2]); 
%
cmap = colormap(hsv(3))
cmap(2,:) = [1 1 1]

subplot(1,2,1); contourf(PX*PY'); title('Truth'); colormap(cmap)
subplot(1,2,2); contourf(R,'LineWidth',0); title('Correlation Matrix'); colormap(cmap)


%
[A B D] = chunkSMR(R,k,'als',Dx,Dy); Rhat = A*diag(D)*B';
[A1 B1 D1] = chunkSMR(R,k,'als'); Rhat1 = A1*diag(D1)*B1';
[u s v] = svds(R,k); 
Rhatsvd = u*s*v';

s = pinv(PX)*R*pinv(PY)'
Rhatoracle = PX*s*PY'; 
plot(Rhatoracle,R); shg
%%
close all; 
% check FIT 
SStot = trace(R'*R); 
SStott = trace(Rt'*Rt); 
E = R - Roracle; SSEoracle = trace(E'*E);
E = Rt - Roracle; SSEoraclet = trace(E'*E);
E = R - Rhat; SSEc = trace(E'*E);
E = Rt - Rhat; SSEt = trace(E'*E);
E = R - Rhatsvd; SSEcsvd = trace(E'*E);
E = Rt - Rhatsvd; SSEtsvd = trace(E'*E);
E = R - Rhat1; SSEc1 = trace(E'*E);
E = Rt - Rhat1; SSEt1 = trace(E'*E);
SSE = [SSEc SSEt; SSEc1 SSEt1; SSEcsvd SSEtsvd;SStot SStott];
cLab = {'Calibrattion';'Validation'};
rLab = {'D corr';'D = 1';'SVD'}
bar(SSE'); set(gca,'Xticklabel',cLab)
legend(rLab); shg

a = A; b = B; rhat = Rhat;
%% 
close all;
cmap = colormap(hsv(3))
cmap = redbluecmap(7)
% cmap(2,:) = [1 1 1]
%
figure
 
subplot(2,2,1); zlm = [-1 1];contourf(PX*PY','LineWidth',0); title('Truth'); colormap(cmap); caxis(zlm)
subplot(2,2,2); zlm = [-1 1]*max(abs(R(:))); contourf(R,'LineWidth',0); title('Correlation Matrix');caxis(zlm)
subplot(2,2,3); zlm = [-1 1]*max(abs(Rhat(:)));contourf(Rhat,'LineWidth',0); title('Estimated cluster Matrix D = corr');caxis(zlm)
subplot(2,2,4); zlm = [-1 1]*max(abs(Rhat1(:))); contourf(Rhat1,'LineWidth',0); title('Estimated cluster Matrix D = FL');caxis(zlm)


figure;
subplot(1,2,1);plot(a*diag(sign(sum(a)))); hold on; plot(colNorm(  PX),':','linewidth',2); hold off; view(-90,90) % swap the x and y axis
subplot(1,2,2);plot(b*diag(sign(sum(b)))); hold on; plot(colNorm(  PY),':','linewidth',2); hold off;
shg
%%
close all; 
zlm = [-1 1];contourf(PX*PY','LineWidth',0); colormap(cmap); caxis(zlm)
% print('paper/Fig5A','-depsc'); close all; 
matlab2tikz('paper/Fig5A.tex'); close all; 
zlm = [-1 1]*max(abs(R(:))); contourf(R,'LineWidth',0); colormap(cmap); caxis(zlm)
% print('paper/Fig5B','-depsc'); close all; 
matlab2tikz('paper/Fig5B.tex'); close all; 
zlm = [-1 1]*max(abs(Rhat(:)));contourf(Rhat,'LineWidth',0); colormap(cmap); caxis(zlm)
% print('paper/Fig5C','-depsc'); close all; 
matlab2tikz('paper/Fig5C.tex'); close all; 
zlm = [-1 1]*max(abs(Rhat1(:))); contourf(Rhat1,'LineWidth',0); colormap(cmap); caxis(zlm)
% print('paper/Fig5D','-depsc'); close all; 
matlab2tikz('paper/Fig5D.tex'); close all; 

% matrix2latex(SSE, 'paper/chunckTable.tex', 'rowLabels', rLab, 'columnLabels', cLab, 'alignment', 'c', 'format', '%-6.1f', 'size', 'tiny'); 




%%
for i=1:k
    P = zeros(p,1);
    id = sort(rand(p,1));
    id = id > sp(i) & id < sp(i+1);
    P(id) = 1*i;
    PP1(:,i) = P;
    P2 = zeros(n,1);
    id = sort(rand(n,1));
    id = id > sp(i) & id < sp(i+1);
    P2(id ) = 1*i;
    PP2(:,i) = P2;
end

% X = P2*P' + randn(n,p)*0.3;
X = PP2*PP1' + randn(n,p)*0.3;
