%% PLOT F1 counts vs F1 theta
colors = {'DarkMagenta','DarkTurquoise'};
for d = 1:2
    
    figure(20);subplot(2,2,3)
    feat1 = POPv3.lick.thetaF1{d};
    feat2 = POPv3.lick.countsF1{d};
    
    hold on;scatter(feat1(:,1),feat2(:,1),'filled','markerfacecolor',rgb(colors{d}));
    set(gca,'xlim',[0 1],'ylim',[0 1],'ytick',0:.5:1,'xtick',0:.5:1);
    xlabel('F1 lick theta');ylabel('F1 lick count');
    hold on; plot([0 1],[0 1],'-.k')
    axis square
    
    figure(20);subplot(2,2,4)
    scatter(feat1(:,2),feat2(:,2),'filled','markerfacecolor',rgb(colors{d}));
    set(gca,'xlim',[0 1],'ylim',[0 1],'ytick',0:.5:1,'xtick',0:.5:1);
    xlabel('F1 nolick theta');ylabel('F1 nolick count');
    hold on; plot([0 1],[0 1],'-.k')
    axis square
end
%% PLOT RMSE IN BOXPLOT FORM
% figure(30);clf;boxplot(cell2mat(POP.RMSE))

figure(30);clf;boxplot([POPv2.logRMSE{2:3}])
set(gca,'xticklabel',{'Counts','Theta','All'})
hold on; plot([3.5; 3.5],[0 0; 0.5 .5],'-.k')
set(gca,'ylim',[0 .5],'ytick',[0:.1:1])
ylabel('RMSE')


figure(650);clf
for i = 1:length(POP.RMSE{1})
    hold on;scatter([1 2 3],POP.RMSE{1}(i,:),'k')
    hold on;plot([1 2 3],POP.RMSE{1}(i,:),'-k')
end
%% PLOTTING RMSE in SCATTER FORM
colors = {'DarkMagenta','DarkTurquoise'};
contvals = POPv2.logRMSEDROP{3};
contvals = [contvals(:,1) contvals(:,3) contvals(:,2)]
semivals = POPv2.logRMSEDROP{2};
semivals = [semivals(:,1) semivals(:,3) semivals(:,2)]
figure(480);clf

hold on; 
plot(1:3,semivals,'-o','Markeredgecolor',rgb(colors{1}),'color',rgb(colors{1}));
set(gca,'xtick',[],'ytick',0:.25:.5,'ylim',[0 .5])
plot(1:3,contvals,'-o','Markeredgecolor',rgb(colors{2}),'color',rgb(colors{2}));
set(gca,'xtick',[],'ytick',0:.25:.5,'ylim',[0 .5])
axis square
%% PLOTTING ACCURACY in SCATTER FORM
colors = {'DarkMagenta','DarkTurquoise'};
figure(485);clf;
subplot(1,2,2)
plot(1:3,POPv2.performanceDROP.accuracy{2}','-o','Markerfacecolor',rgb(colors{2}),'Markeredgecolor',rgb(colors{2}),'color',rgb(colors{2}));
set(gca,'xtick',[],'ytick',0:25:100,'ylim',[50 100])
axis square
subplot(1,2,1)
plot(1:3,POPv2.performanceDROP.accuracy{1}','-o','Markerfacecolor',rgb(colors{1}),'Markeredgecolor',rgb(colors{1}),'color',rgb(colors{1}));
set(gca,'xtick',[],'ytick',0:25:100,'ylim',[50 100])
axis square
%% PLOTTING F1 in SCATTER FORM
colors = {'DarkMagenta','DarkTurquoise'};
lickcont = [POPv2.lickDROP.countsF1{2}(:,1) POPv2.lickDROP.thetaF1{2}(:,1) POPv2.lickDROP.ubered{2}(:,1)];
nolickcont = [POPv2.lickDROP.countsF1{2}(:,2) POPv2.lickDROP.thetaF1{2}(:,2) POPv2.lickDROP.ubered{2}(:,2)];
licksemi = [POPv2.lickDROP.countsF1{1}(:,1) POPv2.lickDROP.thetaF1{1}(:,1) POPv2.lickDROP.ubered{1}(:,1)];
nolicksemi = [POPv2.lickDROP.countsF1{1}(:,2) POPv2.lickDROP.thetaF1{1}(:,2) POPv2.lickDROP.ubered{1}(:,2)];

figure(481);clf
subplot(2,2,2)
plot(1:3,lickcont,'-o','Markerfacecolor',rgb(colors{2}),'Markeredgecolor',rgb(colors{2}),'color',rgb(colors{2}));
set(gca,'xtick',[],'ytick',0:.25:1,'ylim',[0 1])
axis square
title('lick')
subplot(2,2,1)
plot(1:3,licksemi,'-o','Markerfacecolor',rgb(colors{1}),'Markeredgecolor',rgb(colors{1}),'color',rgb(colors{1}));
set(gca,'xtick',[],'ytick',0:.25:1,'ylim',[0 1])
title('lick')
axis square
subplot(2,2,4)
plot(1:3,nolickcont,'-o','Markerfacecolor',rgb(colors{2}),'Markeredgecolor',rgb(colors{2}),'color',rgb(colors{2}));
set(gca,'xtick',[],'ytick',0:.25:1,'ylim',[0 1])
title('withhold lick')
axis square
subplot(2,2,3)
plot(1:3,nolicksemi,'-o','Markerfacecolor',rgb(colors{1}),'Markeredgecolor',rgb(colors{1}),'color',rgb(colors{1}));
set(gca,'xtick',[],'ytick',0:.25:1,'ylim',[0 1])
title('withhold lick')
axis square
%% PLOTTING RAW FEATURES in scatter form LOG
feat = POPv2.logfeatures.ubered;
figure(470);clf
hold on;
for k = 1:length(feat)
    scatter(feat{k}(1,:)./feat{k}(3,:),feat{k}(2,:)./feat{k}(3,:),'filled','markerfacecolor',rgb(colors{k}))
end
set(gca,'xtick',[-5:1:0],'xlim',[-2 0])
xlabel('raw theta weight bias normalized');ylabel('raw counts weight bias normalized')

%% F1 SCORE DEVIATION (actual F1 - ideal F1)
colors = {'DarkMagenta','DarkTurquoise'};
alltheta = cell2mat(POPv2.lick.thetaF1)-cell2mat(POPv2.gng.thetaF1);
allcounts = cell2mat(POPv2.lick.countsF1)-cell2mat(POPv2.gng.countsF1);

droppedtheta = cell2mat(POPv2.lickDROP.thetaF1)-cell2mat(POPv2.gngDROP.thetaF1);
droppedcounts = cell2mat(POPv2.lickDROP.countsF1)-cell2mat(POPv2.gngDROP.countsF1);


idealthetadiff = cell2mat(POPv2.gngDROP.thetaF1)-cell2mat(POPv2.gng.thetaF1);
idealcountsdiff = cell2mat(POPv2.gngDROP.countsF1)-cell2mat(POPv2.gng.countsF1);
actualthetadiff = cell2mat(POPv2.lickDROP.thetaF1)-cell2mat(POPv2.lick.thetaF1);
actualcountsdiff = cell2mat(POPv2.lickDROP.countsF1)-cell2mat(POPv2.lick.countsF1);

thetatoplot = idealthetadiff;
countstoplot = idealcountsdiff;

figure(39);subplot(2,1,1)
scatter(thetatoplot(:,1),countstoplot(:,1),'filled','markerfacecolor',rgb(colors{1}))
hold on; 
scatter(thetatoplot(:,3),countstoplot(:,3),'filled','markerfacecolor',rgb(colors{2}))
plot([0 0],[-1 1],'-.k');
plot([-1 1],[0 0],'-.k');
set(gca,'xlim',[-1 1],'ylim',[-1 1])
axis square

subplot(2,1,2)
scatter(thetatoplot(:,2),countstoplot(:,2),'filled','markerfacecolor',rgb(colors{1}))
hold on; 
scatter(thetatoplot(:,4),countstoplot(:,4),'filled','markerfacecolor',rgb(colors{2}))
plot([0 0],[-1 1],'-.k');
plot([-1 1],[0 0],'-.k');
set(gca,'xlim',[-1 1],'ylim',[-1 1])
axis square

%% F1 SCORE DEVIATION (actual F1 - ideal F1) x SBIAS
colors = {'DarkMagenta','DarkTurquoise'};
allF1modidx = (cell2mat(POPv2.lick.thetaF1) - cell2mat(POPv2.lick.countsF1)) ./ (cell2mat(POPv2.lick.thetaF1) + cell2mat(POPv2.lick.countsF1));
touchF1modidx = (cell2mat(POPv2.lickDROP.thetaF1) - cell2mat(POPv2.lickDROP.countsF1)) ./ (cell2mat(POPv2.lickDROP.thetaF1) + cell2mat(POPv2.lickDROP.countsF1));


diffvals = touchF1modidx-allF1modidx;



figure(28040);clf
subplot(2,2,1)
scatter(POPv2.SBIAS{2},allF1modidx(:,1),'filled','markerfacecolor',rgb(colors{1}))
hold on; scatter(POPv2.SBIAS{3},allF1modidx(:,3),'filled','markerfacecolor',rgb(colors{2}))
set(gca,'xlim',[-1 0],'ylim',[-1 1])
subplot(2,2,3)
scatter(POPv2.SBIAS{2},allF1modidx(:,2),'filled','markerfacecolor',rgb(colors{1}))
hold on; scatter(POPv2.SBIAS{3},allF1modidx(:,4),'filled','markerfacecolor',rgb(colors{2}))
set(gca,'xlim',[-1 0],'ylim',[-1 1])
subplot(2,2,2)
scatter(POPv2.SBIAS{2},touchF1modidx(:,1),'filled','markerfacecolor',rgb(colors{1}))
hold on; scatter(POPv2.SBIAS{3},touchF1modidx(:,3),'filled','markerfacecolor',rgb(colors{2}))
set(gca,'xlim',[-1 0],'ylim',[-1 1])
subplot(2,2,4)
scatter(POPv2.SBIAS{2},touchF1modidx(:,2),'filled','markerfacecolor',rgb(colors{1}))
hold on; scatter(POPv2.SBIAS{3},touchF1modidx(:,4),'filled','markerfacecolor',rgb(colors{2}))
set(gca,'xlim',[-1 0],'ylim',[-1 1])



figure(28043);clf
subplot(2,1,1)
scatter(POPv2.SBIAS{2},diffvals(:,1),'filled','markerfacecolor',rgb(colors{1}))
hold on; scatter(POPv2.SBIAS{3},diffvals(:,3),'filled','markerfacecolor',rgb(colors{2}))
set(gca,'xlim',[-1 0],'ylim',[-1 1])
subplot(2,1,2)
scatter(POPv2.SBIAS{2},diffvals(:,2),'filled','markerfacecolor',rgb(colors{1}))
hold on; scatter(POPv2.SBIAS{3},diffvals(:,4),'filled','markerfacecolor',rgb(colors{2}))
set(gca,'xlim',[-1 0],'ylim',[-1 1])
%% SBIAS vs TREE FEATURE IMPORTANCE
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
popfeatdom = [];
figure(38);clf
for d = 2:3
    
    xax = POPv2.SBIAS{d};
    THETAS = POPv2.logfeatures.ubered{d-1}(1,:)';
    COUNTS = POPv2.logfeatures.ubered{d-1}(2,:)';
    %     THETAS = POP.feature{d}(:,1);
    %     COUNTS = POP.feature{d}(:,2);
    
    featureDom = (abs(THETAS)-abs(COUNTS))./(abs(THETAS)+abs(COUNTS));%multiply by negative 1 so we trend up insteda of down. LOL
    
    figure(38);
    hold on; h = scatter(xax,featureDom,'filled');
    h.CData = rgb(colors{d});
    popfeatdom = [popfeatdom ;featureDom];
    
    
    
end
plot([-1 1],[0 0],'-.k')


% sbias = cell2mat(POP.taskD');

sbias = cell2mat(POPv2.SBIAS);
sbias = [sbias(11:30)]';
ydata = popfeatdom;
[~, orders] = sort(sbias);

[coeff, ~ , mu] = polyfit(sbias(orders),ydata(orders),1);
f = polyval(coeff,sbias(orders),[],mu);

hold on; plot(sbias(orders),f,'k');

modelvals = fitlm(sbias,ydata);
% legend('Discrete','Semi-Continuous','Continuous')
legend('Semi-Continuous','Continuous','location','southeast')
pval = modelvals.Coefficients{2,4};
ordrsq = modelvals.Rsquared.Ordinary;

set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[-1 0],'xtick',-1:.5:1)
%   set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[0 40],'xtick',0:10:40)
%   set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[-1 0],'xtick',-1:1:1) %  EBIAS V3
disp(['Rsquared = ' num2str(ordrsq) ' and pvalue = ' num2str(pval)])
xlabel('Exploration bias')
ylabel('Counts more predictive ------ Theta more predictive')
set(figure(38), 'Units', 'pixels', 'Position', [0, 0, 600, 750]);


%% change in weights
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
popfeatdom = [];
vals=[];
figure(38);clf
figure(380);clf
for d = 2:3
    
    xax = POPv3.SBIAS{d};
    THETAS = POPv3.logfeatures.ubered{d-1}(1,:)';
    COUNTS = POPv3.logfeatures.ubered{d-1}(2,:)';
    
%     dropTHETAS = POPv3.logfeaturesDROP.ubered{d-1}(1,:)';
%     dropCOUNTS = POPv3.logfeaturesDROP.ubered{d-1}(2,:)';
    
    
    featureDom = (abs(THETAS)-abs(COUNTS))./(abs(THETAS)+abs(COUNTS));%multiply by negative 1 so we trend up insteda of down. LOL
    featureDomdrop = (abs(dropTHETAS)-abs(dropCOUNTS))./(abs(dropTHETAS)+abs(dropCOUNTS));
    
    figure(38); hold on;
    plot(1:2,[featureDom featureDomdrop],'-o','Markerfacecolor',rgb(colors{d}),'Markeredgecolor',rgb(colors{d}),'color',rgb(colors{d}))
    set(gca,'xtick',[],'ytick',[-1 0 1],'ylim',[-1 1],'xlim',[.90 2.1])
    vals = [vals ; featureDomdrop-featureDom];
    
    figure(380);hold on
    scatter(xax,featureDomdrop-featureDom,'filled','markerfacecolor',rgb(colors{d}))
end

sbias = cell2mat(POPv2.SBIAS);
sbias = [sbias(11:30)]';
ydata = vals;
[~, orders] = sort(sbias);

[coeff, ~ , mu] = polyfit(sbias(orders),ydata(orders),1);
f = polyval(coeff,sbias(orders),[],mu);

figure(380);hold on
plot(sbias(orders),f,'k');
modelvals = fitlm(sbias,ydata);
modelvals.Coefficients(2,4) ;
% legend('Discrete','Semi-Continuous','Continuous')
legend('Semi-Continuous','Continuous','location','northeast')
pval = modelvals.Coefficients{2,4};
ordrsq = modelvals.Rsquared.Ordinary;
disp(['Rsquared = ' num2str(ordrsq) ' and pvalue = ' num2str(pval)])

set(gca,'xlim',[-1 0],'ylim',[0 2],'ytick',[0:1:2],'xtick',-1:.5:1)




