%% Task D x SBIAS CORRELATIOn 


xdataraw = [POP.SBIAS{2:3}]';
% xdata = cell2mat(gngtaskease);

% ydata = cell2mat(gngcountsratios)*-1
% ydata = cell2mat(gngfullthetaranges);
% ydata = cell2mat(POPthetaoverlap');
% ydataraw = cell2mat(POP.lick.countsF1');
ydataraw = cell2mat(POPv2.lick.countsF1');
ydataplot = ydataraw(:,2); 
ydataregress= ydataraw(:,2);
xdataregress=xdataraw;
xdataregress(isnan(ydataregress)) =[];
ydataregress(isnan(ydataregress))=[];


[~, orders] = sort(xdataregress);

[coeff, ~ , mu] = polyfit(xdataregress(orders),ydataregress(orders),1);
g = polyval(coeff,xdataregress(orders),[],mu);

modelvals = fitlm(xdataregress,ydataregress);
pval = modelvals.Coefficients{2,4};
ordrsq = modelvals.Rsquared.Ordinary;
adjrsq = modelvals.Rsquared.Adjusted;
slope = modelvals.Coefficients{2,1};

disp(['Rsquared = ' num2str(ordrsq) ' and pvalue = ' num2str(pval) ' and slope = ' num2str(slope)])



colors = {'DarkMagenta','DarkTurquoise'};
ranges = [1:10;11:20];
figure(23);subplot(2,2,4)
for d = 1:length(colors)
hold on; h=scatter(xdataraw(ranges(d,:)),ydataplot(ranges(d,:)),'filled');
h.CData = rgb(colors{d});
end
hold on; plot(xdataregress(orders),g,'k')

set(gca,'ylim',[0 1],'ytick',-1:.5:1,'xlim',[-1 0],'xtick',[-1:1:1])
xlabel('Exploration bias') 
set(figure(23), 'Units', 'pixels', 'Position', [0, 0, 600, 600]);
% legend('Semi-Continuous','Continuous','location','southwest')