%% build touch onset curves 
clear Q

for rec=1:length(U)
    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    touchIdx = touchIdx(touchIdx<(numel(spikes)-151)); %elim last touches
    spikesAligned = zeros(numel(touchIdx),201);
    
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-50:150]);
    end
    
    touchONN{rec}=spikesAligned;
    Q.touchraw{rec}=spikesAligned;
end


%% identify window after touch onset from when spiking returns to baseline

cellmean=cellfun(@mean, touchONN,'uniformoutput',0);
baselineFR=zeros(length(cellmean),1);
stdFR=baselineFR;
APstart=stdFR; %GREEN time post touch that AP rate reaches max ap rate/x(set at 3)
b2b=stdFR; %RED time post peak touch that AP rate returns back to baseline+sd
HRpeak=stdFR;%BLUE: time from post peak touch for AP rate to return to peak/x(set to 3 right now)

for d = 1:length(baselineFR)
    
    baselineFRtmp=mean(cellmean{d}(25:50));
    stdFR(d)=std(cellmean{d}(25:50));
    baselineFR(d)=baselineFRtmp+stdFR(d); %baselineFR+std of FR
    
    Rpeak=50+find(cellmean{d}(51:150)==max(cellmean{d}(51:150)),1); %max response peak  after touch onset
    APstart(d)=50+find(cellmean{d}(51:Rpeak)>(max(cellmean{d}(51:150))/3),1);
    HRpeak(d)=Rpeak+find(cellmean{d}(Rpeak:150)<(max(cellmean{d}(51:150))/3),1);
    b2b(d)=Rpeak+(find(cellmean{d}(Rpeak:end)<=baselineFR(d),1)); %time point from -50 of when FR returns back to baseline
        
end


%% OPTIONAL plot window to look at spikes window
figure(39);clf;
plotrow=6;
plotcol=7;

for k = 1:length(touchONN)
    subplot(plotrow,plotcol,k);
    bar(-50:150,sum(touchONN{k})/size(touchONN{k},1),'k');
    hold on; plot([b2b(k)-51 b2b(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'r','linewidth',2)
    hold on; plot([HRpeak(k)-51 HRpeak(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'linewidth',2)
    hold on; plot([APstart(k)-51 APstart(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'g','linewidth',2)
    
    set(gca,'xlim',[-50 150]);
    if k==1
        xlabel('Time from all touch onset (ms)')
        ylabel('spks / ms')
    end
end

%% build PSTH around touch index based on time from touch to back to baseline
values = cell(1,length(U));
grp=values;
grpraw=values;

for rec = 1:length(U)
    wndowSTART=APstart(rec)-51; %choose window to start looking at spikes, if want to look at from when touch starts, use 50
    wndowEND=HRpeak(rec)-51; %choose window to stop looking at spikes 

    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    
    % aside to only choose go touches 
    touchnID=ceil(touchIdx/4000);  
    goT=find(U{rec}.meta.trialType==1);
    [tf, loc] = ismember(touchnID, goT);%find touches that are in go trials
    %touchIdx=touchIdx(tf); %only look at touches in go trials so can easily compare between cont and semi-cont
   
    for var = [1]
        theta = squeeze(U{rec}.S_ctk(var,:,:));
        thetaAtTouch = theta(touchIdx);
        thetaAtTouch = thetaAtTouch(~isnan(thetaAtTouch));
        touchIdx = touchIdx(~isnan(thetaAtTouch));
        spkswndow=spikes(repmat(touchIdx,1,wndowEND-wndowSTART+1)+repmat([wndowSTART:wndowEND],numel(thetaAtTouch),1));%build matrix of spikes post T depending on time to return back to baseline post touch
        spksraw=spikes(repmat(touchIdx,1,51)+repmat([0:50],numel(thetaAtTouch),1)); %raw spikes  built around touch... need to incorporate still 6/2/17
        ISI=cell(1,size(spkswndow,1));
        for i = 1:size(spkswndow,1)
            spks=find(spkswndow(i,:)==1);
            spksdif=spks-[0 spks(1:end-1)];
            ISI{i}=spksdif(2:end);
        end
        ISI=cellfun(@nanmean,ISI);
        spikesAtTouch = sum(spkswndow,2);
        baselinePreTouch = sum(spikes(repmat(touchIdx,1,25)+repmat([-25:-1],numel(touchIdx),1)),2);
        grp{rec}=[thetaAtTouch spikesAtTouch baselinePreTouch ISI'];
        grpraw{rec}=[thetaAtTouch spksraw];
    end
end

%% Touch Quant - mod IDX and ISI 
% Builds a modulation index that measures change in AP rate from 25ms pre
% touch to b2b val post touch
%final product modIdx
% col 1 = theta 
% col 2 = spk rate post touch (determined by b2b)
% col 3 = spk rate pre touch 
% col 4 = avg ISI post touch
% col 5 = avg spk count 
% col 6 = modulation idx (AP rate post - AP rate pre)/(AP rate post + AP
% rate pre).
Q.valNames = {'theta', 'AP rate post touch', 'AP rate pre touch', 'ISI post touch', 'touch count', 'mod idx'};
Q.TW = [APstart'-51; HRpeak'-51; b2b'-51];
Q.baseline = baselineFR;
modIdx=cell(length(U),1);
close all;

for rec=1:length(U)
    
    [Q.theta.sorted{rec}, Q.theta.sortedBy{rec}, Q.theta.binBounds{rec}]=binslin(grp{rec}(:,1),grp{rec},'equalE',41,-50,50);%1 degree bins
    touchspks=cell2mat(cellfun(@(x) nanmean(x,1),Q.theta.sorted{rec},'uniformoutput',0));
    touchspks(:,5)=cellfun(@numel, Q.theta.sorted{rec});
    selected=find(touchspks(:,5)>nansum(touchspks(:,5))*.015); %selecting only for thetas that have at least 1.5% of total touches.
    touchspks(:,1)=linspace(mean(Q.theta.binBounds{rec}(1:2)),mean(Q.theta.binBounds{rec}(end-1:end)),length(Q.theta.binBounds{rec})-1);
    
    [Q.thetaraw.sorted{rec}, Q.thetaraw.sortedBy{rec}, Q.thetaraw.binBounds{rec}]=binslin(grpraw{rec}(:,1),grpraw{rec},'equalE',41,-50,50);%1 degree bins
    
    histocoll=cell2mat(cellfun(@(x) nanmean(x,1),Q.thetaraw.sorted{rec},'uniformoutput',0));
    histocoll(:,1)=linspace(mean(Q.thetaraw.binBounds{rec}(1:2)),mean(Q.thetaraw.binBounds{rec}(end-1:end)),length(Q.thetaraw.binBounds{rec})-1);
    
    Q.rawbinned{rec}=histocoll; %raw spikes from 0:50ms after touch
    Q.selected{rec}=selected; %selected rows based on hvaing at least 1.5% of total touches
    
    modIdxEQ = (touchspks(:,2)-touchspks(:,3))./(touchspks(:,2)+touchspks(:,3));%spks post touch - spks 25ms pre touch/ spks post touch + spks25mspretouch
    modIdxEQ(isnan(modIdxEQ))=0;
    modIdx{rec}=[touchspks modIdxEQ];
    
    pmodidx=(max(modIdxEQ(selected))-min(modIdxEQ(selected)))/(max(modIdxEQ(selected))+min(modIdxEQ(selected)));
    
    figure(10);subplot(6,7,rec)
    plot(modIdx{rec}(selected,6),modIdx{rec}(selected,1))
    %text(.7,modIdx{rec}(selected(end-3),1),num2str(pmodidx));
    set(gca,'ylim',[min(modIdx{rec}(selected,1)) max(modIdx{rec}(selected,1))],'xlim',[-1 1]);
    
    figure(12);subplot(6,7,rec)
    %plot(modIdx{rec}(selected,2),modIdx{rec}(selected,1))
     plot(modIdx{rec}(selected,4),modIdx{rec}(selected,1),'r');
    %text(.7,modIdx{rec}(selected(end-3),1),num2str(pmodidx));
    set(gca,'ylim',[min(modIdx{rec}(selected,1)) max(modIdx{rec}(selected,1))])
    
    
%     figure(11);subplot(6,7,rec)
%     plot(modIdx{rec}(selected,4),modIdx{rec}(selected,1));
%     set(gca,'ylim',[min(modIdx{rec}(selected,1)) max(modIdx{rec}(selected,1))]);
%     
end
Q.val = modIdx;
%% Touch cell or naw

for d=1:length(U)

postFR=mean(mean(Q.touchraw{d}(:,[Q.TW(1,d)+50:Q.TW(2,d)+50])));
%postFR=mean(mean(Q.touchraw{d}(:,[5+50:20+50])));
tmod(d)=(postFR-Q.baseline(d))/(postFR+Q.baseline(d));
tcell(d)=postFR>Q.baseline(d)*1.25;
end

ton = find(abs(tmod)>.1);

figure(39);clf;
plotrow=6;
plotcol=7;

for k = 1:length(touchONN)
    subplot(plotrow,plotcol,k);
    bar(-50:150,sum(touchONN{k})/size(touchONN{k},1),'k');
    %hold on; plot([b2b(k)-51 b2b(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'r','linewidth',2)
    %hold on; plot([HRpeak(k)-51 HRpeak(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'linewidth',2)
    %hold on; plot([APstart(k)-51 APstart(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'g','linewidth',2)
    
    if intersect(k,ton)>0
        text(75,max(sum(touchONN{k})/size(touchONN{k},1)),'******','color','red')
    end
    set(gca,'xlim',[-50 150]);
    if k==1
        xlabel('Time from all touch onset (ms)')
        ylabel('spks / ms')
    end
end
%% 4 moments
for rec = 1:length(U)
    [sorted, sortedBy, binBounds]=binslin(grpraw{rec}(:,1),grpraw{rec},'equalE',41,-50,50);%1 degree bins
    
    histocoll=cell2mat(cellfun(@(x) nanmean(x,1),sorted,'uniformoutput',0));
    histocoll(:,1)=linspace(mean(binBounds(1:2)),mean(binBounds(end-1:end)),length(binBounds)-1);
    
end










