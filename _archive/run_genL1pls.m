cd '~/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/MATLAB/work/Sparsity_in_chemometrics/genSMR'
clear; clc; close all;
load '~/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/Course and teaching/Advanced Chemometrics 2013-2014/data/beer.mat'

% lambda1 = 1; 
% lambda2 = 1; 
% D = [flD*lambda2; lD*lambda1]; 
x = beer.data;
mx = mean(x);
sx = sqrt(std(x));
x = (x -  ones(size(x,1),1)*mx)*diag(1./sx); 
%
%%x = scale(x,mx,sx); 
% xt = scale(beertest.data,mx,sx);
xt = (beertest.data -  ones(size(beertest.data,1),1)*mx)*diag(1./sx); 
my = mean(extract.data);
sy = std(extract.data);
y = (extract.data - my)/sy;
yt = (extracttest.data - my)/sy;
% [y my sy] = auto(extract.data); 
% yt = scale(extracttest.data,my,sy);
wl = beer.axisscale{2};
save NIRbeer.mat x xt y yt wl my sy mx sx
%%
% [T P W B] = genL1pls(x,y,2,D);
% crossvalidate it 
results = cvgenL1pls(x,y,3,[150 300],20,1);
results2 = cvgenL1pls(x,y,3,[150 300],20,0);

close all; 
[regP idbest] = plot_cvgenL1pls(results2);

nlv = 5; 
resSIMpls = crossval(x,y,'sim',results.cvid,12);
[reg,ssq,xlds,ylds,wts,xscrs,yscrs,basis] = simpls(x,y,nlv);

res = results
clc; 
[regP idbest] = plot_cvgenL1pls(res);
[regP2 idbest2] = plot_cvgenL1pls(results2);
res.design(:,idbest)
res.designlb
close all;
xax = beer.axisscale{2};
b = [reg(nlv,:); regP';regP2']*10; 
save genPLS_beerResults.mat
%%
load genPLS_beerResults.mat
figure; 
plot(sum(abs(results.B)>1e-7),results.RMSE,'o'); ylim([0 2]); shg
plot(results.design(3,:), sum(abs(results.B)>1e-7),'o'); shg

subplot(1,2,1);
plot(results.B(end,:),  results.RMSE,'o'); text(results.B(end,:),  results.RMSE,num2str((1:length(results.RMSE))'))
subplot(1,2,2);
plot(results.B(1,:),  results.RMSE,'o'); text(results.B(end,:),  results.RMSE,num2str((1:length(results.RMSE))'))
%%
clear; clc;
close all; 
legFS = 14;
load genPLS_beerResults.mat
% Figure 1
plot(xax,b); 
col = get(gca,'colororder')
hold on; 
plot(xax,beer.data','color',ones(1,3)*0.5); 
h = legend(['PLS (' num2str(nlv) 'LV)'],'L1 constrained PLS (3LV)(*)','L1 constrained PLS (3LV)','location','southeast');
set(h,'FontSize',legFS)
set(gca,'ytick',0)
hold off; 
xlabel('wavelenght (nm)')
shg; 
% matlab2tikz('paper/Fig3A.tex'); close all; 
print('paper/Fig3A_v2','-depsc'); close all; 

figure; 
yobs = y*sy + my;
yobstest = yt*sy + my; 
% PLS
leg = []; 
plot(yobs,resSIMpls.cvpred(:,1,nlv)*sy + my,'o','color',col(1,:)); 
leg = [leg; {'PLS cv'}];
e = yobs - (resSIMpls.cvpred(:,1,nlv)*sy + my); rmsecv_pls = sqrt(e'*e / length(e));
hold on; 
plot(yobstest,sy*xt*reg(nlv,:)' + my,'*','color',col(1,:)); 
leg = [leg; {'PLS testset'}];
e = yobstest - (sy*xt*reg(nlv,:)'+ my); rmsep_pls = sqrt(e'*e / length(e));
% penPLS LS updated
plot(yobs,res.Yhat(:,idbest)*sy + my,'o','color',col(2,:)); 
leg = [leg; {'PLS constr. LS ud cv'}];
e = yobs - (res.Yhat(:,idbest)*sy + my); rmsecv_penpls = sqrt(e'*e / length(e));
plot(yobstest,sy*xt*regP + my,'*','color',col(2,:)); 
leg = [leg; {'PLS constr. LS ud testset'}];
e = yobstest - (sy*xt*regP + my); rmsep_penpls = sqrt(e'*e / length(e));

% penPLS 
plot(yobs,results2.Yhat(:,idbest2)*sy + my,'o','color',col(3,:)); 
leg = [leg; {'PLS constr. ud cv'}];
e = yobs - (results2.Yhat(:,idbest2)*sy + my); rmsecv_penpls2 = sqrt(e'*e / length(e));
plot(yobstest,sy*xt*regP2 + my,'*','color',col(3,:)); 
leg = [leg; {'PLS constr. ud testset'}];
e = yobstest - (sy*xt*regP2 + my); rmsep_penpls2 = sqrt(e'*e / length(e));

hold off; 
abline(1,0,'color','k');
xlabel('measured'); ylabel('predicted')
h2 = legend(leg,'location','southeast'); 
set(h2,'FontSize',legFS)

%matlab2tikz('paper/Fig3B.tex'); close all; 
print('paper/Fig3B_v2','-depsc'); close all; 
%%
{'CV'; 'Test'}
[rmsecv_pls rmsep_pls; rmsecv_penpls rmsep_penpls; rmsecv_penpls2 rmsep_penpls2]
{'PLS'; 'genPLS LSupdate'; 'genPLS'}
