
figure(30);clf;boxplot(cell2mat(POP.RMSE))
set(gca,'xticklabel',{'Counts','Theta','All'})
hold on; plot([3.5 6.5; 3.5 6.5],[0 0; 0.5 .5],'-.k')
set(gca,'ylim',[0 .4],'ytick',[0:.1:.4])
ylabel('RMSE')

%% SBIAS vs TREE FEATURE IMPORTANCE
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
popfeatdom = [];
figure(38);clf
for d = 1:3
    
    xax = POP.SBIAS{d};
    THETAS = POP.feature{d}(:,1);
    COUNTS = POP.feature{d}(:,2);
    
    featureDom = (abs(THETAS)-abs(COUNTS))./(abs(THETAS)+abs(COUNTS));
    
    figure(38);
    hold on; h = scatter(xax,featureDom,'filled');
    h.CData = rgb(colors{d});
    popfeatdom = [popfeatdom ;featureDom];
   
    
    
end
plot([-40 40],[0 0],'-.k')


% sbias = cell2mat(POP.taskD');
% ydata = popfeatdom;

sbias = cell2mat(POP.SBIAS');
ydata = popfeatdom;
[~, orders] = sort(sbias);

[coeff, ~ , mu] = polyfit(sbias(orders),ydata(orders),1);
f = polyval(coeff,sbias(orders),[],mu);
hold on; plot(sbias(orders),f,'k');

modelvals = fitlm(sbias,ydata);
modelvals.Coefficients(2,4) ;
legend('Discrete','Semi-Continuous','Continuous')
pval = modelvals.Coefficients{2,4};
adjrsq = modelvals.Rsquared.Adjusted;

 set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[-20 20])
disp(['Rsquared = ' num2str(adjrsq) ' and pvalue = ' num2str(pval)])
xlabel('Search Bias')
ylabel('Counts More Predictive ------ Theta More Predictive')
set(figure(38), 'Units', 'pixels', 'Position', [0, 0, 600, 750]);