function learningcurvesPlotter(mousel)
%% Plotting Raw Behavioral Data - This will provide the curves for each individual mouse.

%DEFAULT PLOTTING PARAMS
smoothWindow = 200; 
perfthresh=.75;
stretch = 1000; %trials to plot post learning threshold

%-----------------------------
%Plotting for all trials 
catmat = nan(10,30000);
for i = 1:length(mousel)
    
    smoothed = [smooth(mousel{i},smoothWindow) ;nan(stretch,1)];
    AccThresh = find(smoothed>perfthresh);
    crossIdx = find(AccThresh>1000,1);
    
    
    learnedtrial = AccThresh(crossIdx);
    backwards = fliplr(smoothed(1:learnedtrial+stretch)');
    
    catmat(i,1:length(backwards)) = backwards;
    learnedtrialall(i) = learnedtrial;
end
forwards = fliplr(catmat);
figure(24);subplot(3,1,[1 2])
plot(1:length(catmat),forwards,'color',[.8 .8 .8])
xlabel('Trials to expert (thousands)')
ylabel('Accuracy (%)')
set(gca,'xlim',[27000 30000],'xtick',[0:500:30000],'xticklabel',fliplr((0-(stretch/1000):.5:30-(stretch/1000))),'ytick',[0:.25:1],'yticklabel',[0:25:100],'ylim',[.25 1])
hold on; plot([0 30000],[perfthresh perfthresh],'-.k') %plotting expert line threshold
hold on; plot([0 30000],[.5 .5],'-.k') %plotting chance line threshold
hold on; plot(1:length(catmat),fliplr(nanmean(catmat)),'k','linewidth',4)


%-----------------------------
%Plotting number of trials required to reach learning threshold
figure(24);subplot(3,1,[3])
allends = learnedtrialall;
allends(2) = []; %removing one outlier of 23000 trials
scatter(allends,ones(1,length(allends)),150,'filled','MarkerFaceColor',[.8 .8 .8]);
hold on; errorbar(mean(allends),1,std(allends),'horizontal','ko','markersize',15,'markerfacecolor','k')
set(gca,'xlim',[0 10000],'xtick',0:2000:10000,'xticklabel',0:2:10,'yticklabel',[],'ytick',[],'xdir','reverse')
xlabel('Trials to expert (thousands)')

