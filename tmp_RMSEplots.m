
%plot RMSE side plots 
tmp = RMSEgroup';
figure(30);bar(tmp')

legend('Counts','Theta','Uber')
camroll(-90)
ylabel('RMSE')

%%
%build and keep RMSE and optimal feature from uber/counts/theta pred lick
KEEP = RMSEgroup;
KEEP2 = optfeat;







%% run uber_maxexcursion to get search bias/mouse
%build plots with 
THETASwt = KEEP2(2,:)'; %from uber and weights
COUNTSwt = KEEP2(3,:)';
SBIAS = distmed';

THETASrmse = KEEP(:,2); %from uber and psycho
COUNTSrmse = KEEP(:,1);

figure(281);clf;subplot(1,3,1);h=scatter(SBIAS,THETASrmse,'filled');
h.CData = rgb('turquoise'); 
hold on;scatter(SBIAS,COUNTSrmse,'filled','b');
xlabel('Search Bias');ylabel('RMSE')
legend('Theta','Counts')
set(gca,'xlim',[-round(max(abs(SBIAS)),-1) round(max(abs(SBIAS)),-1)]) 

subplot(1,3,2);h=scatter(SBIAS,abs(THETASwt),'filled');
h.CData = rgb('turquoise'); 
hold on;scatter(SBIAS,abs(COUNTSwt),'filled','b');
xlabel('Search Bias');ylabel('Absolute Weight')
legend('Theta','Counts')
title([U{rec}.meta.layer ' Global Comparison']) 
set(gca,'xlim',[-round(max(abs(SBIAS)),-1) round(max(abs(SBIAS)),-1)]) 

subplot(1,3,3);scatter(SBIAS,abs(THETASwt)./abs(COUNTSwt),'filled','r')
ylabel('Theta:Counts Weight');
xlabel('Search Bias')
set(gca,'xlim',[-round(max(abs(SBIAS)),-1) round(max(abs(SBIAS)),-1)]) 

set(figure(281), 'Units', 'pixels', 'Position', [0, 0, 1700, 800]);
print(figure(281),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_globalCompare' ])

%%
counts = RMSEdecomp{1};
thetas = RMSEdecomp{2};
ubered = RMSEdecomp{3};

countsBetter = (counts-thetas)<0;%negative values means theta is larger RMSE and so counts better
thetaBetter = (counts-thetas)>0; %positive vals means theta is smaller RMSE and thus better val


figure (3); clf
for rec = 1:size(counts,2)
    real=[U{rec}.meta.motorPosition;V(rec).trialNums.matrix(5,:)]';%only taking lick row and motorPos row
    [realsorted]= binslin(real(:,1),real,'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
    reallickmean=cell2mat(cellfun(@(x) mean(x,1),realsorted,'uniformoutput',0));
    reallickstd=cell2mat(cellfun(@(x) std(x,0,1),realsorted,'uniformoutput',0));
    
    bestModel = [counts(:,rec),thetas(:,rec)];
    [~ ,x] = min(bestModel,[],2);
    
    figure(3);subplot(2,5,rec)
    plot(xranges',reallickmean(:,2),'r','linewidth',2)
    hold on; scatter(xranges(x==1)',reallickmean(x==1,2),'filled','b');
    hold on; h = scatter(xranges(x==2)',reallickmean(x==2,2),'filled');
    h.CData = rgb('turquoise');
%     hold on; h = scatter(xranges(x==3)',reallickmean(x==3,2),'filled');
%     h.CData = rgb('goldenrod');
    
     set(gca,'xtick',[xranges(1) xranges(6) xranges(11)],'xticklabel',[-1 0 1],'ylim',[0 1],'xlim',[xranges(1)-5000 xranges(end)+5000])
        
    if rec == 5
        legend ('Mouse Psychometric Curve','Counts Better','Theta Better','Uber Better','location','southeast')
    end
end













