function [bestB id] = plot_cvgenL1pls(results)

figure; 
lv = max(results.design(1,:));
idplot = 2; 
unlmb = unique(results.design(idplot,:)); 
nr = floor(sqrt(length(unlmb))); 
nc = ceil(length(unlmb) / nr); 
ylm = [0 max(results.RMSE)];
for i=1:length(unlmb); 
    subplot(nr,nc,i)
    ic = results.design(idplot,:)==unlmb(i);
    d = results.design(:,ic);
    rmse = reshape(results.RMSE(ic),lv,sum(ic)/lv); 
    d = reshape(d(3,:),lv,sum(ic)/lv); 
    plot(d(1,:),rmse,'-o'); ylim(ylm)
%     set(gca,'xscale','log')
    title([char(results.designlb(idplot)) ' = ' num2str(unlmb(i))])
end

figure
subplot(1,2,1)
[mi id] = min(results.RMSE);
bestB = results.B(:,id); 
plot(results.y,results.Yhat(:,id),'o'); 
try ablinenow; end
try abline(1,0); end
subplot(1,2,2)
plot(results.B(:,id))
