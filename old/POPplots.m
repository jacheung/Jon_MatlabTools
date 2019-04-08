
% figure(30);clf;boxplot(cell2mat(POP.RMSE))

figure(30);clf;boxplot([POP.RMSE{2:3}])
set(gca,'xticklabel',{'Counts','Theta','All'})
hold on; plot([3.5; 3.5],[0 0; 0.5 .5],'-.k')
set(gca,'ylim',[0 .3],'ytick',[0:.1:.4])
ylabel('RMSE')


figure(650);clf
for i = 1:length(POP.RMSE{1})
    hold on;scatter([1 2 3],POP.RMSE{1}(i,:),'k')
hold on;plot([1 2 3],POP.RMSE{1}(i,:),'-k')
end

%% SBIAS vs TREE FEATURE IMPORTANCE
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
popfeatdom = [];
figure(38);clf
for d = 2:3
    
    xax = POP.SBIAS{d};
    THETAS = POP.feature{d}(:,1);
    COUNTS = POP.feature{d}(:,2);
    
    featureDom = (abs(THETAS)-abs(COUNTS))./(abs(THETAS)+abs(COUNTS))*-1;%multiply by negative 1 so we trend up insteda of down. LOL
    
    figure(38);
    hold on; h = scatter(xax,featureDom,'filled');
    h.CData = rgb(colors{d});
    popfeatdom = [popfeatdom ;featureDom];
   
    
    
end
plot([-1 1],[0 0],'-.k')


% sbias = cell2mat(POP.taskD');

sbias = cell2mat(POP.SBIAS)';
sbias = [sbias(11:30)];
ydata = popfeatdom; 
% ydata = ydata(11:30)*-1; 
[~, orders] = sort(sbias);

[coeff, ~ , mu] = polyfit(sbias(orders),ydata(orders),1);
coeff
f = polyval(coeff,sbias(orders),[],mu);
residuals = ydata-f;
SSresid = sum(residuals.^2);
SStotal = length(ydata)-1 * var(ydata);
rs1 = 1-SSresid/SStotal

% hold on; plot(sbias(orders),f,'color',rgb(colors{2}));



hold on; plot(sbias(orders),f,'k');

modelvals = fitlm(sbias,ydata);
modelvals.Coefficients(2,4) ;
% legend('Discrete','Semi-Continuous','Continuous')
legend('Semi-Continuous','Continuous','location','southeast')
pval = modelvals.Coefficients{2,4};
adjrsq = modelvals.Rsquared.Adjusted;

%  set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[0 1],'xtick',0:.5:1)
%   set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[0 40],'xtick',0:10:40)
  set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[-1 1],'xtick',-1:1:1) %  EBIAS V3
disp(['Rsquared = ' num2str(adjrsq) ' and pvalue = ' num2str(pval)])
xlabel('Task ease')
ylabel('Theta more predictive ------ Counts more predictive')
set(figure(38), 'Units', 'pixels', 'Position', [0, 0, 600, 750]);