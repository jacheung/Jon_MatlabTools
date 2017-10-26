


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
    set(gca,'ylim',[-1 1],'ytick',[-1:.5:1])
end
plot([-20 20],[0 0],'-.k')


sbias = cell2mat(POP.SBIAS');
ydata = popfeatdom;

sbias = POP.SBIAS{3};
ydata = featureDom;
[~, orders] = sort(sbias);

[coeff, ~ , mu] = polyfit(sbias(orders),ydata(orders),1);
f = polyval(coeff,sbias(orders),[],mu);
hold on; plot(sbias(orders),f,'k');

modelvals = fitlm(sbias,ydata);
modelvals.Coefficients(2,4) ;
legend('Discrete','Semi-Continuous','Continuous')
pval = modelvals.Coefficients{2,4};
adjrsq = modelvals.Rsquared.Adjusted;

disp(['Rsquared = ' num2str(adjrsq) ' and pvalue = ' num2str(pval)])
xlabel('Search Bias')
ylabel('Counts More Predictive ------ Theta More Predictive')