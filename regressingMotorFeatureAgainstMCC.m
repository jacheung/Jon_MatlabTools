figure(40);clf


for i = 1:15
    
    poleOnset = round(mean(BV{i}.meta.poleOnset*1000))+750;
    thetas = squeeze(BV{i}.S_ctk(1,1:poleOnset,:));
    averages = nanmean(thetas');
    averages = averages-nanmean(averages(1:poleOnset-750));
    sems = nanstd(thetas');
    figure(40);subplot(3,5,i);
    filex_shadedErrorBar(1:length(averages),averages,sems)
    set(gca,'xlim',[0 poleOnset],'ylim',[-20 60])
    hold on;plot([poleOnset-750 poleOnset-750],[-10 10],'-.k')
    
    ppoles = nan(BV{i}.k,50);
    for b = 1:BV{i}.k
        touchIdx = [find(BV{i}.S_ctk(9,:,b)==1) find(BV{i}.S_ctk(12,:,b)==1)];
        if ~isempty(touchIdx)
            tmpThetas = BV{i}.S_ctk(1,touchIdx,b);
            ppoles(b,1:length(tmpThetas)) = tmpThetas;
        end
    end
    % FOR CONTINUOUS/SEMI
    polyinputs = sortrows([BV{i}.meta.motorPosition'  ppoles(:,1)]);
    polyinputs(isnan(polyinputs(:,2)),:)=[];
    [coeff, ~ , mu] = polyfit(polyinputs(:,1),polyinputs(:,2),2);
    
    dbtheta(i) = polyval(coeff,mean(BV{i}.meta.ranges),[],mu);
    restmean(i) = mean(averages(1:poleOnset));
    reststd(i) = std(averages(1:poleOnset));
    
end

[numwWhisks,whiskTimes] = numWhiskPreTouch(BV);
%% ASIDE FOR PRETOUCH WHISK FEATS
% How does whisking change in the preWhisk period

for i = 1:15
    tnums = whiskTimes.touchTrials.trialNumbers{i};
    preTWhiskTimes = whiskTimes.touchTrials.preFirstTouch{i};
    
    mpWT = nan(length(preTWhiskTimes),20);
    ampWT = nan(length(preTWhiskTimes),20);
    for b = 1:length(preTWhiskTimes)
        midpoint = BV{i}.S_ctk(4,:,tnums(b));
        amplitude = BV{i}.S_ctk(3,:,tnums(b));
        preTWhiskTimes{b}(preTWhiskTimes{b}==0)=[];
        if ~isempty(preTWhiskTimes{b})
            mpWT(b,1:length(preTWhiskTimes{b})) = midpoint(preTWhiskTimes{b});
            ampWT(b,1:length(preTWhiskTimes{b})) = amplitude(preTWhiskTimes{b});
        end
        
    end
    
    
    x = repmat(1:20,size(mpWT,1),1);
    rawmp = [x(:) mpWT(:)];
    rawamp = [x(:) ampWT(:)];
    
    rawmp = rawmp(~isnan(rawmp(:,2)),:);
    rawamp =  rawamp(~isnan(rawamp(:,2)),:);
    
    whiskFeat.midpoint.raw{i} = rawmp; 
    whiskFeat.midpoint.mean(i,:) = nanmean(mpWT);
    whiskFeat.midpoint.sem(i,:)  = nanstd(mpWT) ./ sqrt(sum(~isnan(mpWT)));
    
    whiskFeat.amplitude.raw{i} = rawamp; 
    whiskFeat.amplitude.mean(i,:) = nanmean(ampWT);
    whiskFeat.amplitude.sem(i,:)  = nanstd(ampWT) ./ sqrt(sum(~isnan(ampWT)));
    
end

figure(8);clf
figure(9);clf
mprho = zeros(1,15);
amprho = zeros(1,15);
whiskCor = zeros(1,15);
for i = 1:15
    
    if ~isempty(whiskFeat.midpoint.raw{i})
        whiskCor(i) = corr(whiskFeat.midpoint.raw{i}(:,2),whiskFeat.amplitude.raw{i}(:,2));
        [rho,mppval] = corr(whiskFeat.midpoint.raw{i});
        mprho(i) = rho(2);
        [rho, amppval] = corr(whiskFeat.amplitude.raw{i});
        amprho(i) = rho(2);
    end
    
    
    figure(8);
    subplot(5,3,i)
    hold on; filex_shadedErrorBar(1:20,whiskFeat.midpoint.mean(i,:),whiskFeat.midpoint.sem(i,:))
    set(gca,'xlim',[1 8],'ylim',[-10 20])
    
    
    
    figure(9);
    subplot(5,3,i)
    hold on; filex_shadedErrorBar(1:20,whiskFeat.amplitude.mean(i,:),whiskFeat.amplitude.sem(i,:))
    set(gca,'xlim',[1 8],'ylim',[0 30])
    
    
    figure(48);
    subplot(5,3,i)
    plot(1:20,normalize_var(whiskFeat.midpoint.mean(i,:),0,1),'b')
    hold on; plot(1:20,normalize_var(whiskFeat.amplitude.mean(i,:),0,1),'g')
    title(num2str(whiskCor(i)));
end
    
figure(8);
suptitle('midpoint')
figure(9);
suptitle('amplitude')
figure(48);
suptitle('midpoint(b) with amplitude(g)')
     
%% Regressing motor features against MCC amp,mp,anglecounts
meanWhisksPreTouch = cellfun(@mean,numwWhisks.touchTrials.preFirstTouch);
meanWhisksPostTouch = cellfun(@mean,numwWhisks.touchTrials.postFirstTouch);
meanWhisksPostTouchNON = cellfun(@mean,numwWhisks.nontouchTrials.postFirstTouch);

mcc = mccpop.fig9B;

angle_from_db = dbtheta-restmean;
fitlm(angle_from_db,mcc(:,6))

xs = [restmean' reststd' dbtheta' angle_from_db' meanWhisksPreTouch' meanWhisksPostTouch'-meanWhisksPostTouchNON' mprho' amprho' whiskCor'];
ys = mcc(:,5:7);

figure(9);clf
mdl = fitlm(whiskCor(1:11),ys(1:11,3))
mdl.plot
xlabel('amp:mp correlation');ylabel('amplitude MCC')






